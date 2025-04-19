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
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/boundary_parameters.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <cmath>
#include <fstream>
#include <memory>
#include <mpi.h>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
triangulationHandler<dim>::triangulationHandler(
  const userInputParameters<dim> &_user_inputs,
  const MGInfo<dim>              &mg_info)
  : user_inputs(&_user_inputs)
{
  if constexpr (dim == 1)
    {
      triangulation = std::make_shared<Triangulation>(
        dealii::Triangulation<dim>::limit_level_difference_at_vertices);
    }
  else
    {
      triangulation = std::make_shared<Triangulation>(
        MPI_COMM_WORLD,
        dealii::Triangulation<dim>::limit_level_difference_at_vertices);
    }

  has_multigrid = mg_info.has_multigrid();
  min_level     = mg_info.get_mg_min_level();
  max_level     = mg_info.get_mg_max_level();
}

template <int dim>
const typename triangulationHandler<dim>::Triangulation &
triangulationHandler<dim>::get_triangulation() const
{
  Assert(triangulation != nullptr, dealii::ExcNotInitialized());
  return *triangulation;
}

template <int dim>
const std::vector<std::shared_ptr<const dealii::Triangulation<dim>>> &
triangulationHandler<dim>::get_mg_triangulation() const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());
  Assert(!coarsened_triangulations.empty(), dealii::ExcNotInitialized());
  return coarsened_triangulations;
}

template <int dim>
const dealii::Triangulation<dim> &
triangulationHandler<dim>::get_mg_triangulation(unsigned int level) const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());
  Assert(!coarsened_triangulations.empty(), dealii::ExcNotInitialized());
  Assert(coarsened_triangulations.size() >= level,
         dealii::ExcMessage(
           "The coarse triangulation set does not contain that specified level"));
  return *coarsened_triangulations[level];
}

template <int dim>
unsigned int
triangulationHandler<dim>::get_n_global_levels() const
{
  Assert(triangulation != nullptr, dealii::ExcNotInitialized());
  return triangulation->n_global_levels();
}

template <int dim>
unsigned int
triangulationHandler<dim>::get_mg_min_level() const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());
  Assert(!coarsened_triangulations.empty(), dealii::ExcNotInitialized());
  return min_level;
}

template <int dim>
unsigned int
triangulationHandler<dim>::get_mg_max_level() const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());
  Assert(!coarsened_triangulations.empty(), dealii::ExcNotInitialized());
  return max_level;
}

template <int dim>
bool
triangulationHandler<dim>::has_setup_multigrid() const
{
  return has_multigrid;
}

template <int dim>
void
triangulationHandler<dim>::generate_mesh()
{
  // TODO (landinjm): Add more generality in selecting mesh types
  if (user_inputs->spatial_discretization.radius != 0.0)
    {
      // TODO (landinjm): Adding assertion about periodic boundary conditions for spheres
      // Generate a sphere

      // TODO (landinjm): Add assertion that the user cannot specify multiple boundary
      // conditions for a hyper_ball geometry
      dealii::GridGenerator::hyper_ball(*triangulation,
                                        dealii::Point<dim>(),
                                        user_inputs->spatial_discretization.radius);
    }
  else
    {
      // TODO (landinjm): Add assertions about periodic boundary conditions for
      // rectangular domains here. Not sure whether it is better to check for assertions
      // here or when we parse user inputs.

      // Generate rectangle
      dealii::GridGenerator::subdivided_hyper_rectangle(
        *triangulation,
        user_inputs->spatial_discretization.subdivisions,
        dealii::Point<dim>(),
        dealii::Point<dim>(user_inputs->spatial_discretization.size));

      // Mark boundaries. This is done before global refinement to reduce the number of
      // cells we have to loop through.
      mark_boundaries();

      // Mark periodicity
      mark_periodic();
    }

    // TODO (landinjm): This be better as a combination of a parameter flag and debug
    // mode. Output triangulation to vtk if in debug mode
#ifdef DEBUG
  export_triangulation_as_vtk("triangulation");
#endif

  // Global refinement
  triangulation->refine_global(user_inputs->spatial_discretization.global_refinement);

  // Create the triangulations for the coarser levels if we have at least one instance of
  // multigrid for any of the fields
  if (!has_multigrid)
    {
      return;
    }

  Assert(triangulation->n_global_levels() > 1,
         dealii::ExcMessage(
           "Multigrid preconditioners require multilevel triangulations"));

  coarsened_triangulations =
    dealii::MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
      *triangulation);

  // TODO (landinjm): p-multigrid
}

template <int dim>
void
triangulationHandler<dim>::export_triangulation_as_vtk(const std::string &filename) const
{
  const dealii::GridOut grid_out;
  std::ofstream         out(filename + ".vtk");
  grid_out.write_vtk(*triangulation, out);
  conditionalOStreams::pout_base() << "Triangulation written to " << filename << ".vtk\n";
}

template <int dim>
void
triangulationHandler<dim>::mark_boundaries() const
{
  const double tolerance = 1e-12;

  // Loop through the cells
  for (const auto &cell : triangulation->active_cell_iterators())
    {
      // Mark the faces (faces_per_cell = 2*dim)
      for (unsigned int face_number = 0;
           face_number < dealii::GeometryInfo<dim>::faces_per_cell;
           ++face_number)
        {
          // Direction for quad and hex cells
          auto direction = static_cast<unsigned int>(std::floor(face_number / 2));

          // Mark the boundary id for x=0, y=0, z=0
          if (std::fabs(cell->face(face_number)->center()(direction) - 0) < tolerance)
            {
              cell->face(face_number)->set_boundary_id(face_number);
            }
          // Mark the boundary id for x=max, y=max, z=max
          else if (std::fabs(cell->face(face_number)->center()(direction) -
                             (user_inputs->spatial_discretization.size[direction])) <
                   tolerance)
            {
              cell->face(face_number)->set_boundary_id(face_number);
            }
        }
    }
}

template <int dim>
void
triangulationHandler<dim>::mark_periodic()
{
  // Add periodicity in the triangulation where specified in the boundary conditions. Note
  // that if one field is periodic all others should be as well.
  for (const auto &[index, boundary_condition] :
       user_inputs->boundary_parameters.boundary_condition_list)
    {
      for (const auto &[component, condition] : boundary_condition)
        {
          for (const auto &[boundary_id, boundary_type] :
               condition.boundary_condition_map)
            {
              if (boundary_type == boundaryCondition::type::PERIODIC)
                {
                  // Skip boundary ids that are odd since those map to the even faces
                  if (boundary_id % 2 != 0)
                    {
                      continue;
                    }

                  // Create a vector of matched pairs that we fill and enforce upon the
                  // constaints
                  std::vector<dealii::GridTools::PeriodicFacePair<
                    typename Triangulation::cell_iterator>>
                    periodicity_vector;

                  // Determine the direction
                  const auto direction =
                    static_cast<unsigned int>(std::floor(boundary_id / dim));

                  // Collect the matched pairs on the coarsest level of the mesh
                  dealii::GridTools::collect_periodic_faces(*triangulation,
                                                            boundary_id,
                                                            boundary_id + 1,
                                                            direction,
                                                            periodicity_vector);

                  // Set constraints
                  triangulation->add_periodicity(periodicity_vector);
                }
            }
        }
    }
}

INSTANTIATE_UNI_TEMPLATE(triangulationHandler)

PRISMS_PF_END_NAMESPACE
