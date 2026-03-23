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
    , Mc(get_user_inputs().user_constants.get_model_constant_double("Mc"))
    , Mn(get_user_inputs().user_constants.get_model_constant_double("Mn"))
    , Kn(get_user_inputs().user_constants.get_model_constant_double("Kn"))
    , center1(
        get_user_inputs().user_constants.get_model_constant_rank_1_tensor("center1"))
    , center2(
        get_user_inputs().user_constants.get_model_constant_rank_1_tensor("center2"))
    , radius1(get_user_inputs().user_constants.get_model_constant_double("radius1"))
    , radius2(get_user_inputs().user_constants.get_model_constant_double("radius2"))
    , matrix_concentration(get_user_inputs().user_constants.get_model_constant_double(
        "matrix_concentration"))
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    double dist1 = point.distance(center1);
    double dist2 = point.distance(center2);

    if (index == 0)
      {
        scalar_value = matrix_concentration;
        scalar_value += 0.5 * (0.125) * (1.0 - std::tanh((dist1 - radius1) / (1.0)));
        scalar_value += 0.5 * (0.125) * (1.0 - std::tanh((dist2 - radius2) / (1.0)));
      }
    else
      {
        scalar_value += 0.5 * (1.0 - std::tanh((dist1 - radius1) / (1.0)));
        scalar_value += 0.5 * (1.0 - std::tanh((dist2 - radius2) / (1.0)));
      }
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_group_id) const override
  {
    const double dt = sim_timer.get_timestep();
    if (solve_group_id == 0) // explicit
      {
        // --- Getting the values and derivatives of the model variables ---

        const ScalarValue c  = variable_list.template get_value<Scalar, OldOne>(0);
        const ScalarGrad  cx = variable_list.template get_gradient<Scalar, OldOne>(0);

        const ScalarValue n  = variable_list.template get_value<Scalar, OldOne>(1);
        const ScalarGrad  nx = variable_list.template get_gradient<Scalar, OldOne>(1);

        // --- Setting the expressions for the terms in the governing equations ---

        const ScalarValue c2 = c * c;
        const ScalarValue c3 = c * c2;
        const ScalarValue c4 = c2 * c2;

        const ScalarValue n2 = n * n;
        const ScalarValue n3 = n * n2;
        const ScalarValue n4 = n2 * n2;
        const ScalarValue n5 = n * n4;

        // Free energy for each phase and their first and second derivatives
        ScalarValue fa  = (-1.6704 - 4.776 * c + 5.1622 * c2 - 2.7375 * c3 + 1.3687 * c4);
        ScalarValue fac = (-4.776 + 10.3244 * c - 8.2125 * c2 + 5.4748 * c3);
        ScalarValue facc = (10.3244 - 16.425 * c + 16.4244 * c2);
        ScalarValue fb   = (5.0 * c2 - 5.9746 * c - 1.5924);
        ScalarValue fbc  = (10.0 * c - 5.9746);
        ScalarValue fbcc = (10.0);

        // Interpolation function and its derivative
        ScalarValue h  = (10.0 * n3 - 15.0 * n4 + 6.0 * n5);
        ScalarValue hn = (30.0 * n2 - 60.0 * n3 + 30.0 * n4);

        // Residual equations
        ScalarGrad  mux  = (cx * ((1.0 - h) * facc + h * fbcc) + nx * ((fbc - fac) * hn));
        ScalarValue eq_c = c;
        ScalarGrad  eqx_c = ((-Mc * dt) * mux);
        ScalarValue eq_n  = (n - (dt * Mn) * (fb - fa) * hn);
        ScalarGrad  eqx_n = ((-dt * Kn * Mn) * nx);

        // --- Submitting the terms for the governing equations ---

        // Terms for the equation to evolve the concentration
        variable_list.set_value_term(0, eq_c);
        variable_list.set_gradient_term(0, eqx_c);

        // Terms for the equation to evolve the order parameter
        variable_list.set_value_term(1, eq_n);
        variable_list.set_gradient_term(1, eqx_n);
      }
    else if (solve_group_id == 1) // postprocess
      {
        // --- Getting the values and derivatives of the model variables ---

        // c
        const ScalarValue c  = variable_list.template get_value<Scalar, Current>(0);
        const ScalarGrad  cx = variable_list.template get_gradient<Scalar, Current>(0);

        // n
        const ScalarValue n  = variable_list.template get_value<Scalar, Current>(1);
        const ScalarGrad  nx = variable_list.template get_gradient<Scalar, Current>(1);

        // --- Setting the expressions for the terms in the postprocessing expressions
        // ---
        const ScalarValue c2 = c * c;
        const ScalarValue c3 = c * c2;
        const ScalarValue c4 = c2 * c2;

        const ScalarValue n2 = n * n;
        const ScalarValue n3 = n * n2;
        const ScalarValue n4 = n2 * n2;
        const ScalarValue n5 = n * n4;

        // Free energy for each phase and their first and second derivatives
        ScalarValue fa = (-1.6704 - 4.776 * c + 5.1622 * c2 - 2.7375 * c3 + 1.3687 * c4);
        ScalarValue fb = (5.0 * c2 - 5.9746 * c - 1.5924);

        // Interpolation function and its derivative
        ScalarValue h = (10.0 * n3 - 15.0 * n4 + 6.0 * n5);

        // The homogeneous free energy
        ScalarValue f_chem = (1.0 - h) * fa + h * fb;

        // The gradient free energy
        ScalarValue f_grad = 0.5 * Kn * nx * nx;

        // The total free energy
        ScalarValue f_tot = f_chem + f_grad;

        // The magnitude of the gradient of c
        ScalarValue mag_grad_c = cx.norm();

        // --- Submitting the terms for the postprocessing expressions ---

        variable_list.set_value_term(2, f_tot);
        variable_list.set_value_term(3, mag_grad_c);
      }
  }

  number             Mc;
  number             Mn;
  number             Kn;
  dealii::Point<dim> center1;
  dealii::Point<dim> center2;
  number             radius1;
  number             radius2;
  number             matrix_concentration;
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

  constexpr unsigned int dim    = 2;
  constexpr unsigned int degree = 1;

  std::vector<FieldAttributes> fields = {
    FieldAttributes("c", Scalar),
    FieldAttributes("n", Scalar),
    FieldAttributes("f_tot", Scalar),
    FieldAttributes("mag_grad_c", Scalar),

  };

  SolveGroup exp_group(
    0,
    Explicit,
    Primary,
    {0, 1},
    make_dependency_set(fields,
                        {"old_1(c)", "grad(old_1(c))", "old_1(n)", "grad(old_1(n))"}));
  SolveGroup              pp_group(1,
                      Explicit,
                      PostProcess,
                                   {2, 3},
                      make_dependency_set(fields, {"n", "grad(n)", "c", "grad(c)"}));
  std::vector<SolveGroup> solve_groups({exp_group, pp_group});

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
