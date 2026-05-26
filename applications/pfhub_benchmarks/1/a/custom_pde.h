// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int degree, typename number>
class CustomPDE : public PDEOperatorBase<2, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, 2, ScalarValue>;
  using ScalarHess  = dealii::Tensor<2, 2, ScalarValue>;
  using VectorValue = dealii::Tensor<1, 2, ScalarValue>;
  using VectorGrad  = dealii::Tensor<2, 2, ScalarValue>;
  using VectorHess  = dealii::Tensor<3, 2, ScalarValue>;
  using PDEOperatorBase<2, degree, number>::get_user_inputs;
  using PDEOperatorBase<2, degree, number>::get_pf_tools;

  // Intentionally hardcoded because this is a benchmark.
  static constexpr double c_alpha = 0.3;
  static constexpr double c_beta  = 0.7;
  static constexpr double rho     = 5.0;
  static constexpr double kappa   = 2.0;
  static constexpr double M       = 5.0;
  static constexpr double epsilon = 0.01;
  static constexpr double c0      = 0.5;

  /**
   * @brief Constructor.
   */
  CustomPDE(const UserInputParameters<2> &_user_inputs, PhaseFieldTools<2> &_pf_tools)
    : PDEOperatorBase<2, degree, number>(_user_inputs, _pf_tools)
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int     &index,
                        [[maybe_unused]] const unsigned int     &component,
                        [[maybe_unused]] const dealii::Point<2> &point,
                        [[maybe_unused]] number                 &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    using std::cos;
    const double x = point[0];
    const double y = point[1];
    if (index == 0) // redundant
      {
        const double squared_term = cos(0.13 * x) * cos(0.087 * y);
        scalar_value =
          c0 + epsilon * (cos(0.105 * x) * cos(0.11 * y) + squared_term * squared_term +
                          (cos(0.025 * x - 0.15 * y) * cos(0.07 * x - 0.02 * y)));
      }
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<2, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer             &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 0) // c
      {
        ScalarValue c_val   = variable_list.template get_value<Scalar, OldOne>(0);
        ScalarGrad  mu_grad = variable_list.template get_gradient<Scalar, OldOne>(1);

        ScalarGrad c_vec = -M * sim_timer.get_timestep() * -mu_grad;

        variable_list.set_value_term(0, c_val);
        variable_list.set_gradient_term(0, -c_vec);
      }
    else if (solve_block_id == 1) // mu
      {
        constexpr double mid    = 0.5 * (c_alpha + c_beta);
        ScalarValue      c_val  = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad       c_grad = variable_list.template get_gradient<Scalar, Current>(0);

        ScalarValue fcV =
          4.0 * rho * (c_val - mid) * (c_val - c_alpha) * (c_val - c_beta);

        ScalarValue mu_val = fcV;
        ScalarGrad  mu_vec = -kappa * c_grad;

        variable_list.set_value_term(1, mu_val);
        variable_list.set_gradient_term(1, -mu_vec);
      }
    else if (solve_block_id == 2) // pp
      {
        ScalarValue c_val  = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad  c_grad = variable_list.template get_gradient<Scalar, Current>(0);

        ScalarValue f_chemical = rho * (c_val - c_alpha) * (c_val - c_alpha) *
                                 (c_val - c_beta) * (c_val - c_beta);
        ScalarValue f_gradient = 0.5 * kappa * c_grad.norm_square();
        ScalarValue f_tot      = f_chemical + f_gradient;
        variable_list.set_value_term(2, f_tot);
      }
  }
};

PRISMS_PF_END_NAMESPACE
