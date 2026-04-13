// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>

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
    , m_well(1.0)
    , kappa(2.0)
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    constexpr number center[12][3] = {
      {0.1, 0.3,  0},
      {0.8, 0.7,  0},
      {0.5, 0.2,  0},
      {0.4, 0.4,  0},
      {0.3, 0.9,  0},
      {0.8, 0.1,  0},
      {0.9, 0.5,  0},
      {0.0, 0.1,  0},
      {0.1, 0.6,  0},
      {0.5, 0.6,  0},
      {1,   1,    0},
      {0.7, 0.95, 0}
    };
    constexpr number rad[12] = {12, 14, 19, 16, 11, 12, 17, 15, 20, 10, 11, 14};
    number           dist    = 0.0;
    for (unsigned int i = 0; i < 12; i++)
      {
        dist = 0.0;
        for (unsigned int dir = 0; dir < dim; dir++)
          {
            dist += (point[dir] - center[i][dir] * 100.0) *
                    (point[dir] - center[i][dir] * 100.0);
          }
        dist = std::sqrt(dist);

        scalar_value += 0.5 * (1.0 - std::tanh((dist - rad[i]) / 1.5));
      }
    scalar_value = std::min(scalar_value, static_cast<number>(1.0));
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_block_id) const override
  {
    if (solve_block_id == 1) // explicit
      {
        ScalarValue n_val  = variable_list.template get_value<Scalar, OldOne>(0);
        ScalarGrad  n_grad = variable_list.template get_gradient<Scalar, OldOne>(0);
        ScalarValue f_well = 4.0 * n_val * (n_val - 1.0) * (n_val - 0.5);
        ScalarValue eq_n   = n_val - sim_timer.get_timestep() * m_well * f_well;
        ScalarGrad  eqx_n  = -sim_timer.get_timestep() * kappa * m_well * n_grad;

        variable_list.set_value_term(0, eq_n);
        variable_list.set_gradient_term(0, eqx_n);
      }
    else if (solve_block_id == 2) // postprocess
      {
        ScalarValue n_val  = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad  n_grad = variable_list.template get_gradient<Scalar, Current>(0);
        ScalarValue f_tot  = 0.0;
        ScalarValue f_chem =
          n_val * n_val * n_val * n_val - 2.0 * n_val * n_val * n_val + n_val * n_val;
        ScalarValue f_grad = 0.0;
        for (unsigned int i = 0; i < dim; i++)
          {
            f_grad += 0.5 * kappa * n_grad[i] * n_grad[i];
          }
        f_tot = f_chem + f_grad;

        variable_list.set_value_term(1,
                                     std::sqrt(n_grad[0] * n_grad[0] +
                                               n_grad[1] * n_grad[1]));
        variable_list.set_value_term(2, f_tot);
      }
  }

  number m_well;
  number kappa;
};

PRISMS_PF_END_NAMESPACE
