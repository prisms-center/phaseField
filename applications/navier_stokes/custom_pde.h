// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/utilities/utilities.h>

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
    if (index == 0 || index == 1)
      {
        if (component == 0 && boundary_id == 0)
          {
            const auto y           = point[1];
            const auto U_m         = 0.15;
            const auto H           = 0.41;
            vector_component_value = 4.0 * U_m * y * (H - y) / (H * H);
          }
      }
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 0)
      {
        if (sim_timer.get_increment() <= 3)
          {
            const VectorValue u_old = variable_list.template get_value<Vector, OldOne>(0);
            const ScalarGrad  p_hash_grad =
              variable_list.template get_gradient<Scalar, OldOne>(3);
            const ScalarValue tau = sim_timer.get_timestep();

            const VectorValue pressure_term = p_hash_grad;
            const VectorValue timestep_term = -u_old / tau;

            const VectorValue residual = timestep_term + pressure_term;

            variable_list.set_value_term(0, -residual);
          }
        else
          {
            const VectorValue u_old = variable_list.template get_value<Vector, OldOne>(0);
            const VectorValue u_old_2 =
              variable_list.template get_value<Vector, OldTwo>(0);
            const ScalarGrad p_hash_grad =
              variable_list.template get_gradient<Scalar, OldOne>(3);
            const ScalarValue tau = sim_timer.get_timestep();

            const VectorValue pressure_term = p_hash_grad;
            const VectorValue timestep_term = (-2.0 * u_old + 0.5 * u_old_2) / tau;

            const VectorValue residual = timestep_term + pressure_term;

            variable_list.set_value_term(0, -residual);
          }
      }
    else if (solve_block_id == 1)
      {
        if (sim_timer.get_increment() <= 3)
          {
            const ScalarValue u_div =
              variable_list.template get_divergence<Vector, Current>(0);
            const ScalarValue tau = sim_timer.get_timestep();

            variable_list.set_value_term(4, u_div / tau);
          }
        else
          {
            const ScalarValue u_div =
              variable_list.template get_divergence<Vector, Current>(0);
            const ScalarValue tau = sim_timer.get_timestep();

            variable_list.set_value_term(4, 1.5 * u_div / tau);
          }
      }
    else if (solve_block_id == 2)
      {
        if (sim_timer.get_increment() <= 3)
          {
            const VectorValue u = variable_list.template get_value<Vector, Current>(0);
            const ScalarValue p_old = variable_list.template get_value<Scalar, OldOne>(2);
            const ScalarValue phi = variable_list.template get_value<Scalar, Current>(4);
            const VectorValue u_star_old =
              variable_list.template get_value<Vector, OldOne>(1);
            const ScalarValue tau = sim_timer.get_timestep();
            const ScalarValue h   = variable_list.get_element_volume();

            const ScalarValue tau_stabilization =
              stabilization_parameter<dim, degree>(tau, h, u_star_old, nu);

            const ScalarValue p      = p_old + phi;
            const VectorValue u_star = u;
            const ScalarValue p_hash = p;

            variable_list.set_value_term(1, u_star);
            variable_list.set_value_term(2, p);
            variable_list.set_value_term(3, p_hash);
            variable_list.set_value_term(5, tau_stabilization);
          }
        else
          {
            const VectorValue u = variable_list.template get_value<Vector, Current>(0);
            const VectorValue u_old = variable_list.template get_value<Vector, OldOne>(0);
            const ScalarValue p_old = variable_list.template get_value<Scalar, OldOne>(2);
            const ScalarValue phi = variable_list.template get_value<Scalar, Current>(4);
            const ScalarValue phi_old =
              variable_list.template get_value<Scalar, OldOne>(4);
            const VectorValue u_star_old =
              variable_list.template get_value<Vector, OldOne>(1);
            const ScalarValue tau = sim_timer.get_timestep();
            const ScalarValue h   = variable_list.get_element_volume();

            const ScalarValue tau_stabilization =
              stabilization_parameter<dim, degree>(tau, h, u_star_old, nu);

            const ScalarValue p      = p_old + phi;
            const VectorValue u_star = 2.0 * u - u_old;
            const ScalarValue p_hash = p + 4.0 / 3.0 * phi - 1.0 / 3.0 * phi_old;

            variable_list.set_value_term(1, u_star);
            variable_list.set_value_term(2, p);
            variable_list.set_value_term(3, p_hash);
            variable_list.set_value_term(5, tau_stabilization);
          }
      }
  }

  void
  compute_lhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 0)
      {
        if (sim_timer.get_increment() <= 3)
          {
            const VectorValue u     = variable_list.template get_value<Vector, LHS>(0);
            const VectorGrad u_grad = variable_list.template get_gradient<Vector, LHS>(0);
            const VectorValue u_star_old =
              variable_list.template get_value<Vector, OldOne>(1);
            const ScalarValue u_star_div =
              variable_list.template get_divergence<Vector, OldOne>(1);
            const ScalarValue tau = sim_timer.get_timestep();

            const VectorValue advection_term = u_star_old * u_grad + 0.5 * u_star_div * u;
            const VectorValue forcing_term;
            const VectorValue timestep_term = u / tau;

            variable_list.set_value_term(0,
                                         timestep_term + advection_term + forcing_term);
            variable_list.set_gradient_term(0, nu * u_grad);
          }
        else
          {
            const VectorValue u     = variable_list.template get_value<Vector, LHS>(0);
            const VectorGrad u_grad = variable_list.template get_gradient<Vector, LHS>(0);
            const VectorValue u_star_old =
              variable_list.template get_value<Vector, OldOne>(1);
            const ScalarValue u_star_div =
              variable_list.template get_divergence<Vector, OldOne>(1);
            const ScalarValue tau = sim_timer.get_timestep();

            const VectorValue advection_term = u_star_old * u_grad + 0.5 * u_star_div * u;
            const VectorValue forcing_term;
            const VectorValue timestep_term = 1.5 * u / tau;

            variable_list.set_value_term(0,
                                         timestep_term + advection_term + forcing_term);
            variable_list.set_gradient_term(0, nu * u_grad);
          }
      }
    else if (solve_block_id == 1)
      {
        const ScalarGrad phi_grad = variable_list.template get_gradient<Scalar, LHS>(4);

        variable_list.set_gradient_term(4, -phi_grad);
      }
  }

  ScalarValue nu = 1.0 / (100.0 * 200.0);
};

PRISMS_PF_END_NAMESPACE
