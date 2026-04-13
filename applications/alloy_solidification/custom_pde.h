// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>

#include <prismspf/utilities/crystal_symmetry.h>

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
    , c_0(get_user_inputs().user_constants.get_model_constant_double("c_0"))
    , U_0(get_user_inputs().user_constants.get_model_constant_double("U_0"))
    , U_off(get_user_inputs().user_constants.get_model_constant_double("U_off"))
    , y_0(get_user_inputs().user_constants.get_model_constant_double("y_0"))
    , epsilon(get_user_inputs().user_constants.get_model_constant_double("epsilon"))
    , k(get_user_inputs().user_constants.get_model_constant_double("k"))
    , lambda(get_user_inputs().user_constants.get_model_constant_double("lambda"))
    , D_tilde(get_user_inputs().user_constants.get_model_constant_double("D_tilde"))
    , V_tilde(get_user_inputs().user_constants.get_model_constant_double("V_tilde"))
    , l_tilde(get_user_inputs().user_constants.get_model_constant_double("l_tilde"))
    , reg_val(get_user_inputs().user_constants.get_model_constant_double("reg_val"))
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
        // For the concentration field, we just set the initial condition
        // to have uniform undercooling everywhere.
        scalar_value = U_0;
      }
    else if (index == 1)
      {
        // For the order parameter, we just place a small seed. Note that
        // the order parameter ranges from -1 to 1 in this model.

        // Center at the origin with some radius
        const dealii::Point<dim> center;
        const double             radius = 5.0;

        // Compute the distance
        double distance = point.distance(center);

        // Apply tanh
        scalar_value = -std::tanh((distance - radius) / std::numbers::sqrt2);
      }
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_block_id) const override

  {
    const double dt = sim_timer.get_timestep();

    // Explicit U and phi evolution
    if (solve_block_id == 0)
      {
        ScalarValue U        = variable_list.template get_value<Scalar, OldOne>(0);
        ScalarGrad  U_grad   = variable_list.template get_gradient<Scalar, OldOne>(0);
        ScalarValue phi      = variable_list.template get_value<Scalar, OldOne>(1);
        ScalarGrad  phi_grad = variable_list.template get_gradient<Scalar, OldOne>(1);
        ScalarValue xi       = variable_list.template get_value<Scalar, OldOne>(2);

        // Compute the interfacial normal vector for the anisotropy function
        // NOTE: We add a small regularization value to avoid division by zero. As this
        // only becomes problematic when phi_grad is small (the bulk where phi is -1 or 1)
        // it has minimal impact on the interfacial dynamics.
        ScalarGrad normal = phi_grad / (phi_grad.norm() + reg_val);

        // Anisotropy term using the crystal symmetry helper function.
        // NOTE: The Symmetries namespace switches between the expanded polynomials and
        // the actual trigonometric function evaluation. For higher order symmetries,
        // there should exist a crossover point in computational cost between the two
        // methods. We have found that to be around ~5-6; however, your environments
        // performance may differ. Feel free to write out your own polynomail if you would
        // like to.
        ScalarValue a_n = 1.0 + epsilon * Symmetries::cos_theta<4>(normal[0], normal[1]);

        // Coefficient before phi
        ScalarValue tau_phi = (1.0 + (1.0 - k) * U) * a_n * a_n;

        // Coefficient before U
        ScalarValue tau_U = ((1.0 + k) - (1.0 - k) * phi) / 2.0;

        // Antitrapping term
        ScalarGrad j_at;
        j_at = (1.0 / (2.0 * std::numbers::sqrt2)) * (1.0 + (1.0 - k) * U) *
               (xi / tau_phi) * normal;

        // grad_phi and grad_U dot product term
        ScalarValue val_term1 = dt * (1.0 + (1.0 - k) * U) * xi / (2.0 * tau_phi * tau_U);
        ScalarValue val_term2 =
          dt * ((1.0 - k) / 2.0) / (tau_U * tau_U) *
          (phi_grad[0] * (D_tilde * ((1.0 - phi) / 2.0) * U_grad[0] + j_at[0]) +
           phi_grad[1] * (D_tilde * ((1.0 - phi) / 2.0) * U_grad[1] + j_at[1]));

        // Define required equations
        ScalarValue eq_U = (U + val_term1 - val_term2);

        ScalarGrad eqx_U =
          (-1.0 * dt * (D_tilde * ((1.0 - phi) / 2.0) * U_grad + j_at) / tau_U);

        ScalarValue eq_phi = phi + (dt * xi / tau_phi);

        variable_list.set_value_term(0, eq_U);
        variable_list.set_gradient_term(0, eqx_U);

        variable_list.set_value_term(1, eq_phi);
      }
    // Explicit xi
    else if (solve_block_id == 1)
      {
        ScalarValue U        = variable_list.template get_value<Scalar, Current>(0);
        ScalarValue phi      = variable_list.template get_value<Scalar, Current>(1);
        ScalarGrad  phi_grad = variable_list.template get_gradient<Scalar, Current>(1);

        // Compute the interfacial normal vector for the anisotropy function
        ScalarGrad normal = phi_grad / (phi_grad.norm() + reg_val);

        // Anisotropy term and its derivative using the crystal symmetry helper functions
        ScalarValue a_n = 1.0 + epsilon * Symmetries::cos_theta<4>(normal[0], normal[1]);
        ScalarValue d_a_n =
          -4.0 * epsilon * Symmetries::sin_theta<4>(normal[0], normal[1]);

        // dimensionless temperature changes
        ScalarValue y   = variable_list.get_q_point_location()[1];
        ScalarValue t   = sim_timer.get_time();
        ScalarValue tep = (y - y_0 - V_tilde * t) / l_tilde;

        // The anisotropy term that enters in to the equation for xi
        ScalarGrad aniso;
        aniso[0] = a_n * a_n * phi_grad[0] - a_n * d_a_n * phi_grad[1];
        aniso[1] = a_n * a_n * phi_grad[1] + a_n * d_a_n * phi_grad[0];

        // Define the terms in the equations
        ScalarValue eq_xi =
          phi - (phi * phi * phi) -
          (lambda * (1.0 - phi * phi) * (1.0 - phi * phi) * (U + tep + U_off));

        ScalarGrad eqx_xi = -aniso;

        variable_list.set_value_term(2, eq_xi);
        variable_list.set_gradient_term(2, eqx_xi);
      }
    // Postprocess c
    else if (solve_block_id == 2)
      {
        ScalarValue U   = variable_list.template get_value<Scalar, Current>(0);
        ScalarValue phi = variable_list.template get_value<Scalar, Current>(1);

        ScalarValue c = (c_0 / 2.0 / (1 + U_0 - U_0 * k)) *
                        ((1.0 + k) - (1.0 - k) * phi) * ((1.0) + (1.0 - k) * U);

        variable_list.set_value_term(3, c);
      }
  }

  number c_0;
  number U_0;
  number U_off;
  number y_0;
  number epsilon;
  number k;
  number lambda;
  number D_tilde;
  number V_tilde;
  number l_tilde;
  number reg_val;
};

PRISMS_PF_END_NAMESPACE
