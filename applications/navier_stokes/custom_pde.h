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
        if (component == 0 && boundary_id == 2)
          {
            const auto y           = point[1];
            const auto U_m         = 1.5;
            const auto H           = 4.1;
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
        const VectorValue u_old = variable_list.template get_value<Vector, OldOne>(0);
        const VectorGrad  u_grad_old =
          variable_list.template get_gradient<Vector, OldOne>(0);
        const ScalarValue u_div_old =
          variable_list.template get_divergence<Vector, OldOne>(0);
        const ScalarValue tau = sim_timer.get_timestep();
        const ScalarValue h   = variable_list.get_element_volume();

        const ScalarValue tau_stabilization =
          stabilization_parameter<dim, degree>(tau, h, u_old, nu);

        const VectorValue advection_term = u_old * u_grad_old + 0.5 * u_div_old * u_old;
        const VectorValue forcing_term;
        const VectorValue timestep_term = -u_old / tau;

        const VectorValue residual_star = timestep_term + advection_term + forcing_term;

        VectorGrad supg_term;
        for (unsigned int i = 0; i < dim; i++)
          {
            for (unsigned int j = 0; j < dim; j++)
              {
                supg_term[i][j] = tau_stabilization * residual_star[i] * u_old[j];
              }
          }

        variable_list.set_value_term(1, -timestep_term - advection_term - forcing_term);
        variable_list.set_gradient_term(1, -supg_term);
      }
    else if (solve_block_id == 1)
      {
        const VectorValue u_old = variable_list.template get_value<Vector, OldOne>(0);
        const VectorGrad  u_grad_old =
          variable_list.template get_gradient<Vector, OldOne>(0);
        const ScalarValue u_div_old =
          variable_list.template get_divergence<Vector, OldOne>(0);
        const VectorValue u_star = variable_list.template get_value<Vector, Current>(1);
        const ScalarValue u_star_div =
          variable_list.template get_divergence<Vector, Current>(1);
        const VectorValue u_star_lap =
          variable_list.template get_laplacian<Vector, Current>(1);
        const ScalarValue tau = sim_timer.get_timestep();
        const ScalarValue h   = variable_list.get_element_volume();

        const ScalarValue tau_stabilization =
          stabilization_parameter<dim, degree>(tau, h, u_old, nu);

        const VectorValue advection_term = u_old * u_grad_old + 0.5 * u_div_old * u_old;
        const VectorValue diffusion_term = -nu * u_star_lap;
        const VectorValue forcing_term;
        const VectorValue timestep_term = (u_star - u_old) / tau;

        const VectorValue residual_star =
          timestep_term + advection_term + diffusion_term + forcing_term;

        variable_list.set_value_term(2, -u_star_div / tau);
        variable_list.set_gradient_term(2, -tau_stabilization * residual_star / tau);
      }
    else if (solve_block_id == 2)
      {
        const VectorValue u_old  = variable_list.template get_value<Vector, OldOne>(0);
        const VectorValue u_star = variable_list.template get_value<Vector, Current>(1);
        const ScalarGrad p_grad = variable_list.template get_gradient<Scalar, Current>(2);
        const ScalarValue tau   = sim_timer.get_timestep();
        const ScalarValue h     = variable_list.get_element_volume();

        const ScalarValue tau_stabilization =
          stabilization_parameter<dim, degree>(tau, h, u_old, nu);

        const VectorValue pressure_term = p_grad;
        const VectorValue timestep_term = -u_star / tau;

        const VectorValue residual_hash = timestep_term + pressure_term;

        VectorGrad supg_term;
        for (unsigned int i = 0; i < dim; i++)
          {
            for (unsigned int j = 0; j < dim; j++)
              {
                supg_term[i][j] = tau_stabilization * residual_hash[i] * u_old[j];
              }
          }

        variable_list.set_value_term(0, -timestep_term - pressure_term);
        variable_list.set_gradient_term(0, -supg_term);
      }
  }

  void
  compute_lhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 0)
      {
        const VectorValue u_old  = variable_list.template get_value<Vector, OldOne>(0);
        const VectorValue u_star = variable_list.template get_value<Vector, LHS>(1);
        const VectorGrad  u_star_grad =
          variable_list.template get_gradient<Vector, LHS>(1);
        const VectorValue u_star_lap =
          variable_list.template get_laplacian<Vector, LHS>(1);
        const ScalarValue tau = sim_timer.get_timestep();
        const ScalarValue h   = variable_list.get_element_volume();

        const ScalarValue tau_stabilization =
          stabilization_parameter<dim, degree>(tau, h, u_old, nu);

        const VectorValue diffusion_term = -nu * u_star_lap;
        const VectorValue timestep_term  = u_star / tau;

        const VectorValue residual_star = timestep_term + diffusion_term;

        VectorGrad supg_term;
        for (unsigned int i = 0; i < dim; i++)
          {
            for (unsigned int j = 0; j < dim; j++)
              {
                supg_term[i][j] = tau_stabilization * residual_star[i] * u_old[j];
              }
          }

        variable_list.set_value_term(1, timestep_term);
        variable_list.set_gradient_term(1, nu * u_star_grad + supg_term);
      }
    else if (solve_block_id == 1)
      {
        const ScalarGrad p_grad = variable_list.template get_gradient<Scalar, LHS>(2);

        variable_list.set_gradient_term(2, p_grad);
      }
    else if (solve_block_id == 2)
      {
        const VectorValue u_old = variable_list.template get_value<Vector, OldOne>(0);
        const VectorValue u     = variable_list.template get_value<Vector, LHS>(0);
        const ScalarValue tau   = sim_timer.get_timestep();
        const ScalarValue h     = variable_list.get_element_volume();

        const ScalarValue tau_stabilization =
          stabilization_parameter<dim, degree>(tau, h, u_old, nu);

        const VectorValue timestep_term = u / tau;

        const VectorValue residual_hash = timestep_term;

        VectorGrad supg_term;
        for (unsigned int i = 0; i < dim; i++)
          {
            for (unsigned int j = 0; j < dim; j++)
              {
                supg_term[i][j] = tau_stabilization * residual_hash[i] * u_old[j];
              }
          }

        variable_list.set_value_term(0, timestep_term);
        variable_list.set_gradient_term(0, supg_term);
      }
  }

  ScalarValue nu = 1.0;
};

PRISMS_PF_END_NAMESPACE
