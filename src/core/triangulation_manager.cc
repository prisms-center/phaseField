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

#include <prismspf/user_inputs/spatial_discretization.h>

#include <prismspf/config.h>

#include <fstream>
#include <memory>
#include <mpi.h>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <>
TriangulationManager<1U>::TriangulationManager(bool _has_multigrid)
  : has_multigrid(_has_multigrid)
  , triangulation(Triangulation<1U>::limit_level_difference_at_vertices)
{}

template <unsigned int dim>
TriangulationManager<dim>::TriangulationManager(bool _has_multigrid)
  : has_multigrid(_has_multigrid)
  , triangulation(MPI_COMM_WORLD, Triangulation<dim>::limit_level_difference_at_vertices)
{}

template <unsigned int dim>
const Triangulation<dim> &
TriangulationManager<dim>::get_triangulation() const
{
  return triangulation;
}

template <unsigned int dim>
const std::vector<std::shared_ptr<const Triangulation<dim>>> &
TriangulationManager<dim>::get_mg_triangulation() const
{
  Assert(!coarsened_triangulations.empty(), dealii::ExcNotInitialized());
  return coarsened_triangulations;
}

template <unsigned int dim>
const Triangulation<dim> &
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
TriangulationManager<dim>::generate_mesh(
  const SpatialDiscretization<dim> &discretization_params)
{
  if (discretization_params.type == TriangulationType::Rectangular)
    {
      discretization_params.rectangular_mesh.generate_mesh(triangulation);
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

  export_triangulation_as_vtk("triangulation");

  // Global refinement
  triangulation.refine_global(discretization_params.global_refinement);

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

#include "core/triangulation_manager.inst"

PRISMS_PF_END_NAMESPACE
