// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>

#include <random>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public PDEOperatorBase<dim, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, ScalarValue>;
  using ScalarHess  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;
  using VectorGrad  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorHess  = dealii::Tensor<3, dim, ScalarValue>;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  /**
   * @brief Constructor.
   */
  explicit CustomPDE(const UserInputParameters<dim> &_user_inputs,
                     PhaseFieldTools<dim>           &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , Mc(get_user_inputs().user_constants.get_double("Mc"))
    , Mn(get_user_inputs().user_constants.get_double("Mn"))
    , Kn(get_user_inputs().user_constants.get_double("Kn"))
    , center1(get_user_inputs().user_constants.get_rank_1_tensor("center1"))
    , center2(get_user_inputs().user_constants.get_rank_1_tensor("center2"))
    , radius1(get_user_inputs().user_constants.get_double("radius1"))
    , radius2(get_user_inputs().user_constants.get_double("radius2"))
    , matrix_concentration(
        get_user_inputs().user_constants.get_double("matrix_concentration"))
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
              unsigned int                         solve_block_id) const override
  {
    const double     dt = sim_timer.get_timestep();
    constexpr double a0 = -1.6704;
    constexpr double a1 = -4.776;
    constexpr double a2 = 5.1622;
    constexpr double a3 = -2.7375;
    constexpr double a4 = 1.3687;
    constexpr double b0 = -1.5924;
    constexpr double b1 = -5.9746;
    constexpr double b2 = 5.0;

    if (solve_block_id == 0) // explicit
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

        ScalarValue fa    = a0 + a1 * c + a2 * c2 + a3 * c3 + a4 * c4;
        ScalarValue fa_c  = a1 + 2 * a2 * c + 3 * a3 * c2 + 4 * a4 * c3;
        ScalarValue fa_cc = 2 * a2 + 6 * a3 * c + 12 * a4 * c2;

        ScalarValue fb    = b0 + b1 * c + b2 * c2;
        ScalarValue fb_c  = b1 + 2 * b2 * c;
        ScalarValue fb_cc = 2 * b2;

        // Interpolation function and its derivative
        ScalarValue h  = 10.0 * n3 - 15.0 * n4 + 6.0 * n5;
        ScalarValue hn = 30.0 * n2 - 60.0 * n3 + 30.0 * n4;

        // RHS equations
        ScalarGrad mu_grad =
          -(cx * ((1.0 - h) * fa_cc + h * fb_cc) + nx * ((fb_c - fa_c) * hn));
        ScalarGrad  dcdt_vec = Mc * -mu_grad;
        ScalarValue dndt_val = -Mn * (fb - fa) * hn;
        ScalarGrad  dndt_vec = Kn * Mn * nx;

        // --- Submitting the terms for the governing equations ---

        // Terms for the equation to evolve the concentration
        variable_list.set_value_term(0, c);
        variable_list.set_gradient_term(0, -dt * dcdt_vec);

        // Terms for the equation to evolve the order parameter
        variable_list.set_value_term(1, n + dt * dndt_val);
        variable_list.set_gradient_term(1, -dt * dndt_vec);
      }
    else if (solve_block_id == 1) // postprocess
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
        ScalarValue fa = a0 + a1 * c + a2 * c2 + a3 * c3 + a4 * c4;
        ScalarValue fb = b0 + b1 * c + b2 * c2;

        // Interpolation function and its derivative
        ScalarValue h = 10.0 * n3 - 15.0 * n4 + 6.0 * n5;

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

PRISMS_PF_END_NAMESPACE
