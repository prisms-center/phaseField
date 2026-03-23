// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/dependencies.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/phase_field_tools.h>
#include <prismspf/core/problem.h>
#include <prismspf/core/simulation_timer.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/input_file_reader.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

using namespace prisms;

template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public PDEOperatorBase<dim, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using ScalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using VectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  /**
   * @brief Constructor.
   */
  explicit CustomPDE(const UserInputParameters<dim> &_user_inputs,
                     PhaseFieldTools<dim>           &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , compliance(get_user_inputs().user_constants.get_model_constant_elasticity_tensor(
        "compliance"))
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {}

  void
  set_nonuniform_dirichlet([[maybe_unused]] const unsigned int       &index,
                           [[maybe_unused]] const unsigned int       &boundary_id,
                           [[maybe_unused]] const unsigned int       &component,
                           [[maybe_unused]] const dealii::Point<dim> &point,
                           [[maybe_unused]] number                   &scalar_value,
                           [[maybe_unused]] number &vector_component_value) const override
  {
    scalar_value           = 0.0;
    vector_component_value = 0.0;
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_group_id) const override
  {
    if (solve_group_id == 1) // linear residual
      {
        ScalarValue dist_from_inclusion = variable_list.get_q_point_location().norm();
        ScalarValue inclusion_radius(10.0);

        VectorGrad transformation_strain;
        using std::tanh;
        ScalarValue strain_value =
          0.01 * (0.5 + 0.5 * tanh(10.0 * (dist_from_inclusion - inclusion_radius)));
        for (unsigned int i = 0; i < dim; i++)
          {
            transformation_strain[i][i] = strain_value;
          }
        VectorGrad stress;
        compute_stress<dim, ScalarValue>(compliance, transformation_strain, stress);
        variable_list.set_gradient_term(0, -stress);
      }
  }

  void
  compute_lhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_group_id) const override
  {
    if (solve_group_id == 1) // linear lhs
      {
        VectorGrad ux = variable_list.template get_symmetric_gradient<Vector, LHS>(0);
        VectorGrad stress;
        compute_stress<dim, ScalarValue>(compliance, ux, stress);
        variable_list.set_gradient_term(0, stress);
      }
  }

  constexpr static unsigned int CIJ_tensor_size = (2 * dim) - 1 + (dim / 3);

  dealii::Tensor<2, CIJ_tensor_size, number> compliance;
};

int
main(int argc, char *argv[])
{
  // Initialize MPI
  dealii::Utilities::MPI::MPI_InitFinalize
    mpi_init(argc, argv, dealii::numbers::invalid_unsigned_int);

  // Restrict deal.II console printing
  dealii::deallog.depth_console(0);

  // Parse the command line options (if there are any) to get the name of the input
  // file
  ParseCMDOptions cli_options(argc, argv);

  constexpr unsigned int dim    = 3;
  constexpr unsigned int degree = 1;

  std::vector<FieldAttributes> fields = {FieldAttributes("u", Vector)};

  SolveGroup linear_solve(1,
                          Linear,
                          Uninitialized,
                          {0},
                          {},
                          make_dependency_set(fields, {"grad(lhs(u))"}));

  std::vector<SolveGroup> solve_groups({linear_solve});

  UserInputParameters<dim>       user_inputs(cli_options.get_parameters_filename());
  PhaseFieldTools<dim>           pf_tools;
  CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);
  Problem<dim, degree, double>   problem(fields,
                                       solve_groups,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);
  problem.solve();

  return 0;
}
