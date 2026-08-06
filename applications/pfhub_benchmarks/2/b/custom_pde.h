// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>

#include <numbers>

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
  static constexpr unsigned int p           = 4;
  static constexpr double       c_alpha     = 0.3;
  static constexpr double       c_beta      = 0.7;
  static constexpr double       rho         = std::numbers::sqrt2;
  static constexpr double       kappa_c     = 3.0;
  static constexpr double       kappa_eta   = 3.0;
  static constexpr double       M           = 5.0;
  static constexpr double       omega       = 1.0;
  static constexpr double       alpha       = 5.0;
  static constexpr double       L           = 5.0;
  static constexpr double       epsilon     = 0.05;
  static constexpr double       c0          = 0.5;
  static constexpr double       epsilon_eta = 0.1;
  static constexpr double       psi         = 1.5;

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
    for (unsigned int i = 0; i < p; ++i)
      {
        if (index == i)
          {
            const double squared_term1 =
              (cos((0.046 + 0.001 * i) * x + (0.0405 + 0.001 * i) * y) *
               cos((0.031 + 0.001 * i) * x - (0.004 + 0.001 * i) * y));
            const double squared_term2 =
              cos((0.01 * i) * x - 4.0) * cos((0.007 + 0.01 * i) * y) +
              cos((0.11 + 0.01 * i) * x) * cos((0.11 + 0.01 * i) * y) +
              psi * squared_term1 * squared_term1;
            scalar_value = epsilon_eta * squared_term2 * squared_term2;
            return;
          }
      }
    if (index == p) // redundant
      {
        const double squared_term1 = cos(0.13 * x) * cos(0.087 * y);
        scalar_value =
          c0 + epsilon * (cos(0.105 * x) * cos(0.11 * y) + squared_term1 * squared_term1 +
                          (cos(0.025 * x - 0.15 * y) * cos(0.07 * x - 0.02 * y)));
        return;
      }
  }

  ScalarValue
  f_alpha(const ScalarValue &c) const
  {
    return rho * rho * (c - c_alpha) * (c - c_alpha);
  }

  ScalarValue
  df_alpha_dc(const ScalarValue &c) const
  {
    return 2.0 * rho * rho * (c - c_alpha);
  }

  ScalarValue
  f_beta(const ScalarValue &c) const
  {
    return rho * rho * (c - c_beta) * (c - c_beta);
  }

  ScalarValue
  df_beta_dc(const ScalarValue &c) const
  {
    return 2.0 * rho * rho * (c - c_beta);
  }

  ScalarValue
  h(const dealii::AlignedVector<ScalarValue> &eta) const
  {
    ScalarValue h_val = 0.0;
    for (unsigned int i = 0; i < p; ++i)
      {
        const ScalarValue &eta_i   = eta[i];
        ScalarValue        eta_i_2 = eta_i * eta_i;
        h_val += eta_i_2 * eta_i * (6.0 * eta_i_2 - 15.0 * eta_i + 10.0);
      }
    return h_val;
  }

  ScalarValue
  dh_deta(const ScalarValue &eta_i) const
  {
    return 30.0 * eta_i * eta_i * (1.0 - eta_i) * (1.0 - eta_i);
  }

  ScalarValue
  sum_eta_sq(const dealii::AlignedVector<ScalarValue> &eta) const
  {
    ScalarValue sum = 0.0;
    for (unsigned int i = 0; i < p; ++i)
      {
        const ScalarValue &eta_i = eta[i];
        sum += eta_i * eta_i;
      }
    return sum;
  }

  ScalarValue
  g(const dealii::AlignedVector<ScalarValue> &eta,
    const ScalarValue                        &sum_eta_sq_val) const
  {
    ScalarValue g_val = 0.0;
    for (unsigned int i = 0; i < p; ++i)
      {
        const ScalarValue &eta_i   = eta[i];
        ScalarValue        eta_i_2 = eta_i * eta_i;
        g_val += eta_i_2 * (1.0 - eta_i) * (1.0 - eta_i);
        g_val += 2.0 * alpha * eta_i_2 * (sum_eta_sq_val - eta_i_2);
      }
    return g_val;
  }

  ScalarValue
  dg_deta(const ScalarValue &eta_i, const ScalarValue &sum_eta_sq_val) const
  {
    ScalarValue dgdeta_i = 2.0 * eta_i * (1.0 - 3.0 * eta_i + 2.0 * eta_i * eta_i);
    dgdeta_i += 4.0 * alpha * eta_i * (sum_eta_sq_val - eta_i * eta_i);
    return dgdeta_i;
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<2, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer             &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    const double dt = sim_timer.get_timestep();
    if (solve_block_id == 0) // explicit
      {
        dealii::AlignedVector<ScalarValue> eta_values(p);
        dealii::AlignedVector<ScalarGrad>  eta_grads(p);
        for (unsigned int i = 0; i < p; ++i)
          {
            eta_values[i] = variable_list.template get_value<Scalar, OldOne>(i);
            eta_grads[i]  = variable_list.template get_gradient<Scalar, OldOne>(i);
          }
        ScalarValue c_val   = variable_list.template get_value<Scalar, OldOne>(p);
        ScalarGrad  mu_grad = variable_list.template get_gradient<Scalar, OldOne>(p + 1);

        dealii::AlignedVector<ScalarValue> dF_deta_val(p);
        dealii::AlignedVector<VectorValue> dF_deta_vec(p);

        const ScalarValue sum_eta_sq_val = sum_eta_sq(eta_values);
        for (unsigned int i = 0; i < p; ++i)
          {
            const ScalarValue dh_deta_i = dh_deta(eta_values[i]);
            dF_deta_val[i]              = (f_beta(c_val) - f_alpha(c_val)) * dh_deta_i +
                             omega * dg_deta(eta_values[i], sum_eta_sq_val);

            dF_deta_vec[i] = -kappa_eta * eta_grads[i];
          }

        ScalarGrad dc_dt_vec = -M * -mu_grad;

        for (unsigned int i = 0; i < p; ++i)
          {
            variable_list.set_value_term(i, eta_values[i] + dt * -L * dF_deta_val[i]);
            variable_list.set_gradient_term(i, -dt * -L * dF_deta_vec[i]);
          }

        variable_list.set_value_term(p, c_val);
        variable_list.set_gradient_term(p, -dt * M * dc_dt_vec);
      }
    else if (solve_block_id == 1) // mu
      {
        dealii::AlignedVector<ScalarValue> eta_values(p);
        for (unsigned int i = 0; i < p; ++i)
          {
            eta_values[i] = variable_list.template get_value<Scalar, Current>(i);
          }
        ScalarValue c_val  = variable_list.template get_value<Scalar, Current>(p);
        ScalarGrad  c_grad = variable_list.template get_gradient<Scalar, Current>(p);

        ScalarValue h_val = h(eta_values);
        ScalarValue mu_val =
          df_alpha_dc(c_val) * (1.0 - h_val) + df_beta_dc(c_val) * h_val;
        VectorValue mu_vec = -kappa_c * c_grad;

        variable_list.set_value_term(p + 1, mu_val);
        variable_list.set_gradient_term(p + 1, -mu_vec);
      }
    else if (solve_block_id == 2) // pp
      {
        dealii::AlignedVector<ScalarValue> eta_values(p);
        dealii::AlignedVector<ScalarGrad>  eta_grads(p);
        for (unsigned int i = 0; i < p; ++i)
          {
            eta_values[i] = variable_list.template get_value<Scalar, Current>(i);
            eta_grads[i]  = variable_list.template get_gradient<Scalar, Current>(i);
          }
        ScalarValue c_val  = variable_list.template get_value<Scalar, Current>(p);
        ScalarGrad  c_grad = variable_list.template get_gradient<Scalar, Current>(p);

        ScalarValue h_val = h(eta_values);
        ScalarValue F_val = f_alpha(c_val) * (1.0 - h_val) + f_beta(c_val) * h_val +
                            omega * g(eta_values, sum_eta_sq(eta_values));
        F_val += 0.5 * kappa_c * (c_grad * c_grad);
        for (unsigned int i = 0; i < p; ++i)
          {
            F_val += 0.5 * kappa_eta * (eta_grads[i] * eta_grads[i]);
          }

        variable_list.set_value_term(p + 2, F_val);
      }
  }
};

PRISMS_PF_END_NAMESPACE
