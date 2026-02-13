// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/triangulation_manager.h>

#include <prismspf/config.h>

#include <fstream>
#include <memory>
#include <mpi.h>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <>
TriangulationManager<1U>::TriangulationManager(bool _has_multigrid)
  : has_multigrid(_has_multigrid)
  , triangulation(dealii::Triangulation<1U>::limit_level_difference_at_vertices)
{}

template <unsigned int dim>
TriangulationManager<dim>::TriangulationManager(bool _has_multigrid)
  : has_multigrid(_has_multigrid)
  , triangulation(MPI_COMM_WORLD,
                  dealii::Triangulation<dim>::limit_level_difference_at_vertices)
{}

template <unsigned int dim>
const typename TriangulationManager<dim>::Triangulation &
TriangulationManager<dim>::get_triangulation() const
{
  return triangulation;
}

template <unsigned int dim>
const std::vector<std::shared_ptr<const dealii::Triangulation<dim>>> &
TriangulationManager<dim>::get_mg_triangulation() const
{
  Assert(!coarsened_triangulations.empty(), dealii::ExcNotInitialized());
  return coarsened_triangulations;
}

template <unsigned int dim>
const dealii::Triangulation<dim> &
TriangulationManager<dim>::get_triangulation(unsigned int relative_level) const
{
  Assert(!coarsened_triangulations.empty(), dealii::ExcNotInitialized());
  Assert(coarsened_triangulations.size() >= relative_level,
         dealii::ExcMessage(
           "The coarse triangulation set does not contain that specified level"));
  return *coarsened_triangulations[relative_level];
}

template <unsigned int dim>
void
TriangulationManager<dim>::generate_mesh(const UserInputParameters<dim> &user_inputs)
{
  const SpatialDiscretization<dim> &discretization_params =
    user_inputs.get_spatial_discretization();
  // TODO (landinjm): Add more generality in selecting mesh types
  if (discretization_params.get_radius() != 0.0)
    {
      // TODO (landinjm): Adding assertion about periodic boundary conditions for spheres
      // Generate a sphere

      // TODO (landinjm): Add assertion that the user cannot specify multiple boundary
      // conditions for a hyper_ball geometry
      dealii::GridGenerator::hyper_ball(triangulation,
                                        dealii::Point<dim>(),
                                        discretization_params.get_radius());
    }
  else
    {
      // TODO (landinjm): Add assertions about periodic boundary conditions for
      // Rectangular domains here. Not sure whether it is better to check for assertions
      // here or when we parse user inputs.

      // Generate rectangle
      dealii::GridGenerator::subdivided_hyper_rectangle(
        triangulation,
        discretization_params.get_subdivisions(),
        discretization_params.get_lower_bound(),
        discretization_params.get_upper_bound());

      // Mark boundaries. This is done before global refinement to reduce the number of
      // cells we have to loop through.
      mark_boundaries(discretization_params);

      // Mark periodicity
      mark_periodic(user_inputs);
    }

    // TODO (landinjm): This be better as a combination of a parameter flag and debug
    // mode. Output triangulation to vtk if in debug mode
#ifdef DEBUG
// if(user_inputs ... output triangulation)
#endif
  {
    export_triangulation_as_vtk("triangulation");
  }
  // Global refinement
  triangulation.refine_global(discretization_params.get_global_refinement());

  // Create the triangulations for the coarser levels if we have multigrid for any of the
  // fields
  reinit();
}

template <unsigned int dim>
void
TriangulationManager<dim>::export_triangulation_as_vtk(const std::string &filename) const
{
  const dealii::GridOut grid_out;
  std::ofstream         out(filename + ".vtk");
  grid_out.write_vtk(triangulation, out);
  ConditionalOStreams::pout_base() << "Triangulation written to " << filename << ".vtk\n";
}

// TODO (fractalsbyx): See if this can be made more efficient
template <unsigned int dim>
void
TriangulationManager<dim>::mark_boundaries(
  const SpatialDiscretization<dim> &discretization_params) const
{
  const double tolerance = 1e-12;

  // Loop through the cells
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // Mark the faces (faces_per_cell = 2*dim)
      for (unsigned int face_number = 0;
           face_number < dealii::GeometryInfo<dim>::faces_per_cell;
           ++face_number)
        {
          // Direction for quad and hex cells
          auto direction = face_number / 2;

          // Lower bound and upper bound
          const double lower_bound = discretization_params.get_lower_bound()[direction];
          const double upper_bound = discretization_params.get_upper_bound()[direction];

          // Mark the boundary id for lower and upper bounds
          if (std::fabs(cell->face(face_number)->center()(direction) - lower_bound) <
                tolerance ||
              std::fabs(cell->face(face_number)->center()(direction) - upper_bound) <
                tolerance)
            {
              cell->face(face_number)->set_boundary_id(face_number);
            }
        }
    }
}

template <unsigned int dim>
void
TriangulationManager<dim>::mark_periodic(const UserInputParameters<dim> &user_inputs)
{
  // Create a little set of boundary ids we've already marked
  std::set<unsigned int> periodic_ids;

  // Add periodicity in the triangulation where specified in the boundary conditions. Note
  // that if one field is periodic all others should be as well.
  // TODO (fractalsbyx): Enforce above condition systematically
  for (const auto &[index, boundary_condition] :
       user_inputs.get_boundary_parameters().get_boundary_condition_list())
    {
      for (const auto &[component, condition] : boundary_condition)
        {
          for (const auto &[boundary_id, boundary_type] :
               condition.get_boundary_condition_map())
            {
              if (boundary_type == BoundaryCondition::Type::Periodic)
                {
                  // Skip boundary ids that are odd since those map to the even faces
                  if (boundary_id % 2 != 0)
                    {
                      continue;
                    }

                  // Skip the id if we've already added periodic boundaries to id
                  if (!periodic_ids.insert(boundary_id).second)
                    {
                      continue;
                    }

                  // Create a vector of matched pairs that we fill and enforce upon the
                  // constaints
                  std::vector<dealii::GridTools::PeriodicFacePair<
                    typename Triangulation::cell_iterator>>
                    periodicity_vector;

                  // Determine the direction
                  const auto direction = boundary_id / 2;

                  // Grab the offset vector from one vertex to another
                  dealii::Tensor<1, dim> offset;
                  offset[direction] =
                    user_inputs.get_spatial_discretization().get_size()[direction];

                  // Collect the matched pairs on the coarsest level of the mesh
                  dealii::GridTools::collect_periodic_faces(triangulation,
                                                            boundary_id,
                                                            boundary_id + 1,
                                                            direction,
                                                            periodicity_vector,
                                                            offset);

                  // Set constraints
                  triangulation.add_periodicity(periodicity_vector);
                }
            }
        }
    }
}

// #include "core/triangulation_manager.inst"

PRISMS_PF_END_NAMESPACE
