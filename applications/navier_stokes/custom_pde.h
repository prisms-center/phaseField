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
  {}

  void
  set_dirichlet([[maybe_unused]] const unsigned int       &index,
                [[maybe_unused]] const unsigned int       &boundary_id,
                [[maybe_unused]] const unsigned int       &component,
                [[maybe_unused]] const dealii::Point<dim> &point,
                [[maybe_unused]] number                   &scalar_value,
                [[maybe_unused]] number &vector_component_value) const override
  {
    if (index == 0)
      {
        if (component == 0 && boundary_id == 2)
          {
            const auto y           = point[1];
            const auto U_m         = 1.5;
            const auto H           = 4.1;
            vector_component_value = 4.0 * U_m * y * (H - y) / (H * H);
          }
        else
          {
            vector_component_value = 0.0;
          }
      }
    else if (index == 4)
      {
        scalar_value = 0.0;
      }
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 1)
      {
        const auto u_old     = variable_list.template get_value<Vector, Current>(0);
        const auto u_old_2   = variable_list.template get_value<Vector, OldOne>(0);
        const auto p_old     = variable_list.template get_value<Scalar, Current>(2);
        const auto phi_old   = variable_list.template get_value<Scalar, Current>(4);
        const auto phi_old_2 = variable_list.template get_value<Scalar, OldOne>(4);

        variable_list.set_value_term(1, 2.0 * u_old - u_old_2);
        variable_list.set_value_term(3,
                                     p_old + 4.0 / 3.0 * phi_old - 1.0 / 3.0 * phi_old_2);
      }
    else if (solve_block_id == 2)
      {
        const auto u_old       = variable_list.template get_value<Vector, OldOne>(0);
        const auto u_old_2     = variable_list.template get_value<Vector, OldTwo>(0);
        const auto u_star      = variable_list.template get_value<Vector, OldOne>(1);
        const auto grad_p_star = variable_list.template get_gradient<Scalar, OldOne>(3);
        const ScalarValue tau  = sim_timer.get_timestep();
        const ScalarValue h    = variable_list.get_element_volume();

        const auto timestep_term = (2.0 * u_old - 0.5 * u_old_2) / tau;

        const auto residual_rhs = timestep_term - grad_p_star;

        const auto tau_stabilization =
          stabilization_parameter<dim, degree>(tau, h, u_star, nu);

        VectorGrad supg_term;
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            supg_term[i][j] = tau_stabilization * residual_rhs[i] * u_star[j];

        variable_list.set_value_term(0, residual_rhs);
        variable_list.set_gradient_term(0, supg_term);
      }
    else if (solve_block_id == 3)
      {
        const auto u_old      = variable_list.template get_value<Vector, OldOne>(0);
        const auto u_old_2    = variable_list.template get_value<Vector, OldTwo>(0);
        const auto u          = variable_list.template get_value<Vector, Current>(0);
        const auto grad_u     = variable_list.template get_gradient<Vector, Current>(0);
        const auto div_u      = variable_list.template get_divergence<Vector, Current>(0);
        const auto laplace_u  = variable_list.template get_laplacian<Vector, Current>(0);
        const auto u_star     = variable_list.template get_value<Vector, OldOne>(1);
        const auto div_u_star = variable_list.template get_divergence<Vector, OldOne>(1);
        const auto grad_p_star = variable_list.template get_gradient<Scalar, OldOne>(3);
        const ScalarValue tau  = sim_timer.get_timestep();
        const ScalarValue h    = variable_list.get_element_volume();

        const auto timestep_term  = (3.0 * u - 4.0 * u_old + u_old_2) / (2.0 * tau);
        const auto advection_term = u_star * grad_u + 0.5 * div_u_star * u;

        const auto residual =
          timestep_term + advection_term - nu * laplace_u + grad_p_star;

        const auto tau_stabilization =
          stabilization_parameter<dim, degree>(tau, h, u_star, nu);

        variable_list.set_value_term(4, 1.5 * div_u / tau);
        variable_list.set_gradient_term(4, tau_stabilization * residual);
      }
    else if (solve_block_id == 4)
      {
        const auto p_old = variable_list.template get_value<Scalar, OldOne>(2);
        const auto phi   = variable_list.template get_value<Scalar, Current>(4);

        variable_list.set_value_term(2, p_old + phi);
      }
  }

  void
  compute_lhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 2)
      {
        const auto u          = variable_list.template get_value<Vector, LHS>(0);
        const auto grad_u     = variable_list.template get_gradient<Vector, LHS>(0);
        const auto laplace_u  = variable_list.template get_laplacian<Vector, LHS>(0);
        const auto u_star     = variable_list.template get_value<Vector, OldOne>(1);
        const auto div_u_star = variable_list.template get_divergence<Vector, OldOne>(1);
        const ScalarValue tau = sim_timer.get_timestep();
        const ScalarValue h   = variable_list.get_element_volume();

        const auto timestep_term  = 1.5 * u / tau;
        const auto advection_term = u_star * grad_u + 0.5 * div_u_star * u;

        const auto residual_lhs = timestep_term + advection_term - nu * laplace_u;

        const auto tau_stabilization =
          stabilization_parameter<dim, degree>(tau, h, u_star, nu);

        VectorGrad supg_term;
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            supg_term[i][j] = tau_stabilization * residual_lhs[i] * u_star[j];

        variable_list.set_value_term(0, timestep_term + advection_term);
        variable_list.set_gradient_term(0, nu * grad_u + supg_term);
      }
    else if (solve_block_id == 3)
      {
        const auto grad_phi = variable_list.template get_gradient<Scalar, LHS>(4);

        variable_list.set_gradient_term(4, -grad_phi);
      }
  }

  ScalarValue nu = 0.0015;
};

PRISMS_PF_END_NAMESPACE
