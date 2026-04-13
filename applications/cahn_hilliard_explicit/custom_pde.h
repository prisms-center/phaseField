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
    , McV(get_user_inputs().user_constants.get_model_constant_double("McV"))
    , KcV(get_user_inputs().user_constants.get_model_constant_double("KcV"))
    , WcV(get_user_inputs().user_constants.get_model_constant_double("WcV"))
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    const dealii::Tensor<1, dim> &mesh_size =
      get_user_inputs().spatial_discretization.rectangular_mesh.size;
    if (index == 0)
      {
        double center[12][3] = {
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
        double rad[12] = {12, 14, 19, 16, 11, 12, 17, 15, 20, 10, 11, 14};
        double dist    = 0.0;
        double sdf     = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 12; i++)
          {
            dist = 0.0;
            for (unsigned int dir = 0; dir < dim; dir++)
              {
                double comp_diff = point[dir] - center[i][dir] * mesh_size[dir];
                dist += comp_diff * comp_diff;
              }
            dist = std::sqrt(dist) - rad[i];
            sdf  = std::min(sdf, dist);
          }
        scalar_value += 0.5 * (1.0 - std::tanh(sdf / 1.5));
      }
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 0) // c
      {
        ScalarValue c   = variable_list.template get_value<Scalar, OldOne>(0);
        ScalarGrad  mux = variable_list.template get_gradient<Scalar, OldOne>(1);

        ScalarValue eq_c  = c;
        ScalarGrad  eqx_c = -McV * sim_timer.get_timestep() * mux;

        variable_list.set_value_term(0, eq_c);
        variable_list.set_gradient_term(0, eqx_c);
      }
    else if (solve_block_id == 1) // mu
      {
        ScalarValue c  = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad  cx = variable_list.template get_gradient<Scalar, Current>(0);

        ScalarValue fcV = WcV * c * (c - 1.0) * (c - 0.5);

        ScalarValue eq_mu  = fcV;
        ScalarGrad  eqx_mu = KcV * cx;

        variable_list.set_value_term(1, eq_mu);
        variable_list.set_gradient_term(1, eqx_mu);
      }
    else if (solve_block_id == 2) // pp
      {
        ScalarValue c  = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad  cx = variable_list.template get_gradient<Scalar, Current>(0);

        ScalarValue f_tot  = 0.0;
        ScalarValue f_chem = c * c * c * c - 2.0 * c * c * c + c * c;
        ScalarValue f_grad = 0.0;

        for (unsigned int i = 0; i < dim; i++)
          {
            f_grad += 0.5 * KcV * cx[i] * cx[i];
          }
        f_tot = f_chem + f_grad;
        variable_list.set_value_term(2, f_tot);
      }
  }

  ScalarValue McV;
  ScalarValue KcV;
  ScalarValue WcV;
};

PRISMS_PF_END_NAMESPACE
