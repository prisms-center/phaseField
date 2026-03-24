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
  using ScalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using ScalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using VectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  /**
   * @brief Constructor.
   */
  explicit CustomPDE(const UserInputParameters<dim> &_user_inputs,
                     PhaseFieldTools<dim>           &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , c0(get_user_inputs().user_constants.get_model_constant_double("c0"))
    , U0(get_user_inputs().user_constants.get_model_constant_double("U0"))
    , U_off(get_user_inputs().user_constants.get_model_constant_double("U_off"))
    , y0(get_user_inputs().user_constants.get_model_constant_double("y0"))
    , epsilon(get_user_inputs().user_constants.get_model_constant_double("epsilon"))
    , k(get_user_inputs().user_constants.get_model_constant_double("k"))
    , lambda(get_user_inputs().user_constants.get_model_constant_double("lambda"))
    , Dtilde(get_user_inputs().user_constants.get_model_constant_double("Dtilde"))
    , Vtilde(get_user_inputs().user_constants.get_model_constant_double("Vtilde"))
    , ltilde(get_user_inputs().user_constants.get_model_constant_double("ltilde"))
    , regval(get_user_inputs().user_constants.get_model_constant_double("regval"))
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
        scalar_value = U0;
      }
    else if (index == 1)
      {
        const std::array<double, 3> center = {
          {0.0, 0.0, 0.0}
        };
        const double rad = 5.0;
        // For the order parameter, we just place a small seed. Note that
        // the order parameter ranges from -1 to 1 in this model.
        double dist = 0.0;
        for (unsigned int dir = 0; dir < dim; dir++)
          {
            const number comp_diff = point[dir] - center[dir] * mesh_size[dir];
            dist += comp_diff * comp_diff;
          }
        dist = std::sqrt(dist);

        scalar_value = -std::tanh((dist - rad) / std::sqrt(2));
      }
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_group_id) const override
  {
    const double dt = sim_timer.get_timestep();
    if (solve_group_id == 0) // explicit U and phi
      {
        // The dimensionless solute supersaturation and its derivatives
        ScalarValue U  = variable_list.template get_value<Scalar, OldOne>(0);
        ScalarGrad  Ux = variable_list.template get_gradient<Scalar, OldOne>(0);

        // The order parameter and its derivatives
        ScalarValue phi  = variable_list.template get_value<Scalar, OldOne>(1);
        ScalarGrad  phix = variable_list.template get_gradient<Scalar, OldOne>(1);

        // The auxiliary parameter and its derivatives
        ScalarValue xi = variable_list.template get_value<Scalar, OldOne>(2);

        // --- Setting the expressions for the terms in the governing equations ---

        // Calculation of interface normal vector
        ScalarValue normgradn = std::sqrt(phix.norm_square());
        ScalarGrad  normal    = phix / (normgradn + regval);

        // The cosine of theta
        ScalarValue cth = normal[0];
        // The sine of theta
        ScalarValue sth = normal[1];
        // The cosine of 4 theta
        ScalarValue c4th = (sth * sth * sth * sth) + (cth * cth * cth * cth) -
                           (6.0 * sth * sth * cth * cth);

        // Anisotropic term
        ScalarValue a_n = 1.0 + (epsilon * c4th);

        // coefficient before phi
        ScalarValue tau_phi = (1.0 + (1.0 - k) * U) * a_n * a_n;

        // coefficient before U
        ScalarValue tau_U = (((1.0 + k) / 2.0) - ((1.0 - k) * phi / 2.0));

        // Antitrapping term
        ScalarGrad j_at;
        j_at = (1.0 / (2.0 * std::sqrt(2.0))) * (1.0 + (1.0 - k) * U) * (xi / tau_phi) *
               normal;

        // grad_phi and grad_U dot product term
        ScalarValue val_term1 = dt * (1.0 + (1.0 - k) * U) * xi / (2.0 * tau_phi * tau_U);
        ScalarValue val_term2 =
          dt * ((1.0 - k) / 2.0) / (tau_U * tau_U) *
          (phix[0] * (Dtilde * ((1.0 - phi) / 2.0) * Ux[0] + j_at[0]) +
           phix[1] * (Dtilde * ((1.0 - phi) / 2.0) * Ux[1] + j_at[1]));

        // Define required equations
        ScalarValue eq_U = (U + val_term1 - val_term2);

        ScalarGrad eqx_U =
          (-1.0 * dt * (Dtilde * ((1.0 - phi) / 2.0) * Ux + j_at) / tau_U);

        ScalarValue eq_phi = phi + (dt * xi / tau_phi);

        // --- Submitting the terms for the governing equations ---

        // Terms for the equation to evolve the concentration
        variable_list.set_value_term(0, eq_U);
        variable_list.set_gradient_term(0, eqx_U);

        // Terms for the equation to evolve the order parameter
        variable_list.set_value_term(1, eq_phi);
      }
    else if (solve_group_id == 1) // xi
      {
        // --- Getting the values and derivatives of the model variables ---

        ScalarValue U    = variable_list.template get_value<Scalar, Current>(0);
        ScalarValue phi  = variable_list.template get_value<Scalar, Current>(1);
        ScalarGrad  phix = variable_list.template get_gradient<Scalar, Current>(1);

        // --- Setting the expressions for the terms in the governing equations ---

        // Calculation of interface normal vector
        ScalarValue normgradn = std::sqrt(phix.norm_square());
        ScalarGrad  normal    = phix / (normgradn + regval);

        // The cosine of theta
        ScalarValue cth = normal[0];
        // The sine of theta
        ScalarValue sth = normal[1];

        // The cosine of 4 theta
        ScalarValue c4th = (sth * sth * sth * sth) + (cth * cth * cth * cth) -
                           (6.0 * sth * sth * cth * cth);
        // The sine of 4 theta
        ScalarValue s4th = (4.0 * sth * cth * cth * cth) - (4.0 * sth * sth * sth * cth);

        // Anisotropic term
        ScalarValue a_n = 1.0 + (epsilon * c4th);

        // gradient energy coefficient, its derivative and square
        ScalarValue a_d = -4.0 * epsilon * s4th;

        // dimensionless temperature changes
        ScalarValue y   = variable_list.get_q_point_location()[1]; // The y-component
        ScalarValue t_n = sim_timer.get_time();                    // The time
        ScalarValue tep = ((y - y0 - Vtilde * t_n) / ltilde);

        // The anisotropy term that enters in to the equation for xi
        ScalarGrad aniso;
        aniso[0] = a_n * a_n * phix[0] - a_n * a_d * phix[1];
        aniso[1] = a_n * a_n * phix[1] + a_n * a_d * phix[0];

        // Define the terms in the equations
        ScalarValue eq_xi =
          phi - (phi * phi * phi) -
          (lambda * (1.0 - phi * phi) * (1.0 - phi * phi) * (U + tep + U_off));

        ScalarGrad eqx_xi = -aniso;

        // --- Submitting the terms for the governing equations ---

        variable_list.set_value_term(2, eq_xi);
        variable_list.set_gradient_term(2, eqx_xi);
      }
    else if (solve_group_id == 2) // c postprocess
      {
        ScalarValue U   = variable_list.template get_value<Scalar, Current>(0);
        ScalarValue phi = variable_list.template get_value<Scalar, Current>(1);

        ScalarValue c = (c0 / 2.0 / (1 + U0 - U0 * k)) * ((1.0 + k) - (1.0 - k) * phi) *
                        ((1.0) + (1.0 - k) * U);

        variable_list.set_value_term(3, c);
      }
  }

  number c0;
  number U0;
  number U_off;
  number y0;
  number epsilon;
  number k;
  number lambda;
  number Dtilde;
  number Vtilde;
  number ltilde;
  number regval;
};

PRISMS_PF_END_NAMESPACE