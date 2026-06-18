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
#include <deal.II/grid/grid_tools_geometry.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/triangulation_manager.h>

#include <prismspf/user_inputs/spatial_discretization.h>

#include <fstream>
#include <mpi.h>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <>
TriangulationManager<1U>::TriangulationManager()
  : triangulation(Triangulation<1U>::limit_level_difference_at_vertices)
{}

template <unsigned int dim>
TriangulationManager<dim>::TriangulationManager()
  : triangulation(MPI_COMM_WORLD,
                  Triangulation<dim>::limit_level_difference_at_vertices,
                  Triangulation<dim>::construct_multigrid_hierarchy)
{}

template <unsigned int dim>
void
TriangulationManager<dim>::init_mg()
{
  coarsened_triangulations =
    dealii::MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
      triangulation);
  std::reverse(coarsened_triangulations.begin(), coarsened_triangulations.end());
  // TODO (landinjm): p-multigrid
}

template <unsigned int dim>
void
TriangulationManager<dim>::clear_mg()
{
  coarsened_triangulations.clear();
}

template <unsigned int dim>
bool
TriangulationManager<dim>::has_mg() const
{
  return !coarsened_triangulations.empty();
}

template <unsigned int dim>
unsigned int
TriangulationManager<dim>::num_levels() const
{
  if (!has_mg())
    {
      return 1;
    }
  return coarsened_triangulations.size();
}

template <unsigned int dim>
const Triangulation<dim> &
TriangulationManager<dim>::get_triangulation() const
{
  return triangulation;
}

template <unsigned int dim>
const dealii::Triangulation<dim> &
TriangulationManager<dim>::get_triangulation(unsigned int relative_level) const
{
  if (!has_mg())
    {
      if (relative_level == 0)
        {
          return triangulation;
        }
      AssertThrow(false, dealii::ExcMessage("Invalid relative level for triangulation."));
    }
  return *(coarsened_triangulations[relative_level]);
}

template <unsigned int dim>
const std::vector<dealii::GridTools::PeriodicFacePair<
  typename dealii::Triangulation<dim>::cell_iterator>> &
TriangulationManager<dim>::get_periodic_face_pairs() const
{
  return periodicity_vector;
}

template <unsigned int dim>
double
TriangulationManager<dim>::get_volume() const
{
  return volume;
}

template <unsigned int dim>
void
TriangulationManager<dim>::generate_mesh(
  const SpatialDiscretization<dim> &discretization_params)
{
  if (discretization_params.type == TriangulationType::Rectangular)
    {
      discretization_params.rectangular_mesh.generate_mesh(triangulation);
      discretization_params.rectangular_mesh.collect_periodic_faces(triangulation,
                                                                    periodicity_vector);
    }
  else if (discretization_params.type == TriangulationType::Spherical)
    {
      discretization_params.spherical_mesh.generate_mesh(triangulation);
    }
  else if (discretization_params.type == TriangulationType::Custom)
    {
      AssertThrow(false, FeatureNotImplemented("Custom Triangulation"));
    }
  else
    {
      AssertThrow(false, UnreachableCode("Invalid TriangulationType"));
    }

  // TODO: Once we move to a pattern of manual initialization after default construction,
  // call this separately using the user_inputs output directory.
  export_triangulation_as_vtk("triangulation");

  // Global refinement
  triangulation.refine_global(discretization_params.global_refinement);

  volume = dealii::GridTools::volume(triangulation);
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

#include "core/triangulation_manager.inst"

PRISMS_PF_END_NAMESPACE
