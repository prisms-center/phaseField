// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>
#include <prismspf/core/type_enums.h>

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
  CustomPDE(const UserInputParameters<dim> &_user_inputs, PhaseFieldTools<dim> &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , McV(get_user_inputs().user_constants.get_double("McV"))
    , KcV(get_user_inputs().user_constants.get_double("KcV"))
    , WcV(get_user_inputs().user_constants.get_double("WcV"))
    , ic_type(get_user_inputs().user_constants.get_int("ic_type"))
    , c0(get_user_inputs().user_constants.get_double("c0"))
    , icamplitude(get_user_inputs().user_constants.get_double("icamplitude"))
    , dist(-1.0, 1.0)
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
    if (index == 0) // redundant
      {
        if (ic_type == 0)
          { // Random number generator (Type std::mt19937_64)
            RNGEngine &rng = get_user_inputs().misc_parameters.rng;
            // noise around c0 with amplitude icamplitude
            scalar_value = c0 + icamplitude * dist(rng);
          }
        else if (ic_type == 1)
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
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    const double dt = sim_timer.get_timestep();
    if (solve_block_id == 0) // c
      {
        ScalarValue c_val   = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad  c_grad  = variable_list.template get_gradient<Scalar, Current>(0);
        ScalarValue mu_val  = variable_list.template get_value<Scalar, Current>(1);
        ScalarGrad  mu_grad = variable_list.template get_gradient<Scalar, Current>(1);

        ScalarValue c_old = variable_list.template get_value<Scalar, OldOne>(0);

        ScalarValue r_c_val = (c_old - c_val);
        VectorValue r_c_vec = (dt * McV * mu_grad);

        ScalarValue r_mu_val = mu_val - WcV * c_val * (c_val - 1.0) * (c_val - 0.5);
        VectorValue r_mu_vec = KcV * c_grad;

        variable_list.set_value_term(0, r_c_val);
        variable_list.set_gradient_term(0, -r_c_vec);
        variable_list.set_value_term(1, r_mu_val);
        variable_list.set_gradient_term(1, -r_mu_vec);
      }
    else if (solve_block_id == 2) // pp
      {
        ScalarValue c  = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad  cx = variable_list.template get_gradient<Scalar, Current>(0);

        ScalarValue f_tot  = 0.0;
        ScalarValue f_chem = c * c * c * c - 2.0 * c * c * c + c * c;
        ScalarValue f_grad = 0.5 * KcV * cx.norm_square();
        f_tot              = f_chem + f_grad;
        variable_list.set_value_term(2, f_tot);
      }
  }

  void
  compute_lhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    const double dt = sim_timer.get_timestep();
    if (solve_block_id == 0) // c
      {
        ScalarValue c_val   = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad  c_grad  = variable_list.template get_gradient<Scalar, Current>(0);
        ScalarValue mu_val  = variable_list.template get_value<Scalar, Current>(1);
        ScalarGrad  mu_grad = variable_list.template get_gradient<Scalar, Current>(1);

        ScalarValue delta_c_val  = variable_list.template get_value<Scalar, Change>(0);
        ScalarGrad  delta_c_grad = variable_list.template get_gradient<Scalar, Change>(0);
        ScalarValue delta_mu_val = variable_list.template get_value<Scalar, Change>(1);
        ScalarGrad delta_mu_grad = variable_list.template get_gradient<Scalar, Change>(1);

        ScalarValue j_c_c_val = (-delta_c_val);
        // VectorValue j_c_dc_vec = ;

        // ScalarValue j_c_mu_val = ;
        VectorValue j_c_mu_vec = (dt * McV * delta_mu_grad);

        ScalarValue j_mu_c_val =
          -WcV * (0.5 + 3.0 * (c_val * (c_val - 1.0))) * delta_c_val;
        VectorValue j_mu_c_vec = KcV * delta_c_grad;

        ScalarValue j_mu_mu_val = delta_mu_val;
        // VectorValue j_mu_mu_vec = ;

        variable_list.set_value_term(0, j_c_c_val);
        variable_list.set_gradient_term(0, -j_c_mu_vec);
        variable_list.set_value_term(1, j_mu_c_val + j_mu_mu_val);
        variable_list.set_gradient_term(1, -j_mu_c_vec);
      }
  }

  ScalarValue McV;
  ScalarValue KcV;
  ScalarValue WcV;

  int                                            ic_type;
  number                                         c0;
  number                                         icamplitude;
  mutable std::uniform_real_distribution<number> dist;
};

PRISMS_PF_END_NAMESPACE
