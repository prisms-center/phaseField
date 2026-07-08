// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>

#include <prismspf/core/dependencies.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/problem.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/type_enums.h>

using namespace prismspf;

template <unsigned int dim>
class ChannelWithCylinder : public Mesh<dim>
{
  using Triangulation = typename Mesh<dim>::Triangulation;

  void
  generate_mesh(Triangulation &triangulation) const override
  {
    dealii::GridGenerator::channel_with_cylinder(triangulation, 0.03, 2, 2.0, true);
  }

  void
  mark_boundaries(Triangulation &triangulation) const override
  {
    // deal.II does this for us
    // TODO: Document what the boundary ids are
  }
};

template <unsigned int dim>
class ChannelWithSquare : public Mesh<dim>
{
  using Triangulation = typename Mesh<dim>::Triangulation;

  void
  generate_mesh(Triangulation &triangulation) const override
  {
    dealii::GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::string   filename = "nsbench2.inp";
    std::ifstream file(filename);
    Assert(file, dealii::ExcFileNotOpen(filename));
    grid_in.read_ucd(file);
  }

  void
  mark_boundaries(Triangulation &triangulation) const override
  {
    // The mesh file does this for us
    // TODO: Document what the boundary ids are
  }
};

int
main(int argc, char *argv[])
{
  // Initialize MPI
  prismspf::MPIInitFinalize mpi_init(argc, argv);

  // Parse the command line options (if there are any) to get the name of the input
  // file
  ParseCMDOptions cli_options(argc, argv);

  constexpr unsigned int dim    = 2;
  constexpr unsigned int degree = 2;

  // TODO: Add documentation
  std::vector<FieldAttributes> fields = {FieldAttributes("u", Vector),
                                         FieldAttributes("u_star", Vector),
                                         FieldAttributes("p"),
                                         FieldAttributes("p_hash"),
                                         FieldAttributes("phi")};

  SolveBlock diffusion;
  diffusion.id            = 0;
  diffusion.field_indices = {0};
  diffusion.solve_type    = Linear;
  diffusion.solve_timing  = Initialized;
  diffusion.dependencies_rhs =
    make_dependency_set(fields,
                        {"old_1(u)", "old_2(u)", "old_1(u_star)", "grad(old_1(p_hash))"});
  diffusion.dependencies_lhs = make_dependency_set(
    fields,
    {"lhs(u)", "grad(lhs(u))", "hess(lhs(u))", "old_1(u_star)", "grad(old_1(u_star))"});

  SolveBlock projection;
  projection.id               = 1;
  projection.field_indices    = {4};
  projection.solve_type       = Linear;
  projection.solve_timing     = Uninitialized;
  projection.dependencies_rhs = make_dependency_set(fields,
                                                    {"u",
                                                     "grad(u)",
                                                     "hess(u)",
                                                     "old_1(u)",
                                                     "old_2(u)",
                                                     "old_1(u_star)",
                                                     "grad(old_1(u_star))",
                                                     "grad(old_1(p_hash))"});
  projection.dependencies_lhs = make_dependency_set(fields, {"grad(lhs(phi))"});

  SolveBlock extrapolation;
  extrapolation.id            = 2;
  extrapolation.field_indices = {1, 2, 3};
  extrapolation.solve_type    = Explicit;
  extrapolation.solve_timing  = Uninitialized;
  extrapolation.dependencies_rhs =
    make_dependency_set(fields, {"old_1(p)", "phi", "u", "old_1(u)", "old_1(phi)"});

  std::vector<SolveBlock> solve_blocks({diffusion, projection, extrapolation});

  UserInputParameters<dim> user_inputs(cli_options.get_parameters_filename());
  ChannelWithSquare<dim>   mesh;
  user_inputs.spatial_discretization.custom_mesh = &mesh;
  user_inputs.spatial_discretization.mesh_type   = Custom;

  PhaseFieldTools<dim>           pf_tools;
  CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);
  Problem<dim, degree, double>   problem(fields,
                                       solve_blocks,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);
  problem.solve();

  return 0;
}
