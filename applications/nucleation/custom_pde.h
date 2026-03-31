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

  double c_avg;
  double McV;
  double MnV;
  double KnV;
  double W_barrier;
  double A0;
  double A2;
  double calmin;
  double B0;
  double B2;
  double cbtmin;
  double k1;
  double k2;
  double tau;
  double epsilon;
  double r_nuc;
  double r_freeze;
  double seeding_duration;
  double interface_coeff;

  /**
   * @brief Constructor.
   */
  explicit CustomPDE(const UserInputParameters<dim> &_user_inputs,
                     PhaseFieldTools<dim>           &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , c_avg(get_user_inputs().user_constants.get_model_constant_double("c_avg"))
    , McV(get_user_inputs().user_constants.get_model_constant_double("McV"))
    , MnV(get_user_inputs().user_constants.get_model_constant_double("MnV"))
    , KnV(get_user_inputs().user_constants.get_model_constant_double("KnV"))
    , W_barrier(get_user_inputs().user_constants.get_model_constant_double("W_barrier"))
    , A0(get_user_inputs().user_constants.get_model_constant_double("A0"))
    , A2(get_user_inputs().user_constants.get_model_constant_double("A2"))
    , calmin(get_user_inputs().user_constants.get_model_constant_double("calmin"))
    , B0(get_user_inputs().user_constants.get_model_constant_double("B0"))
    , B2(get_user_inputs().user_constants.get_model_constant_double("B2"))
    , cbtmin(get_user_inputs().user_constants.get_model_constant_double("cbtmin"))
    , k1(get_user_inputs().user_constants.get_model_constant_double("k1"))
    , k2(get_user_inputs().user_constants.get_model_constant_double("k2"))
    , tau(get_user_inputs().user_constants.get_model_constant_double("tau"))
    , epsilon(get_user_inputs().user_constants.get_model_constant_double("epsilon"))
    , r_nuc(get_user_inputs().user_constants.get_model_constant_double("r_nuc"))
    , r_freeze(get_user_inputs().user_constants.get_model_constant_double("r_freeze"))
    , seeding_duration(
        get_user_inputs().user_constants.get_model_constant_double("seeding_duration"))
    , interface_coeff(std::sqrt(2.0 * KnV / W_barrier))
  {}

  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    // Initial condition for the concentration field
    if (index == 0)
      {
        scalar_value = c_avg;
      }
    // Initial condition for the structural order parameter field
    else
      {
        scalar_value = 0.0;
      }
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_block_id) const override
  {
    if (solve_block_id == 0) // explicit
      {
        // --- Getting the values and derivatives of the model variables ---

        // The concentration and its derivatives
        ScalarValue c  = variable_list.template get_value<Scalar, OldOne>(0);
        ScalarGrad  cx = variable_list.template get_gradient<Scalar, OldOne>(0);

        // The order parameter and its derivatives
        ScalarValue n  = variable_list.template get_value<Scalar, OldOne>(1);
        ScalarGrad  nx = variable_list.template get_gradient<Scalar, OldOne>(1);

        // --- Setting the expressions for the terms in the governing equations ---

        // Interpolation function and its derivative
        ScalarValue hV  = 3.0 * n * n - 2.0 * n * n * n;
        ScalarValue hnV = 6.0 * n - 6.0 * n * n;

        double delta_t = sim_timer.get_timestep();

        // KKS model c_alpha and c_beta as a function of c and h
        ScalarValue c_alpha =
          (B2 * (c - cbtmin * hV) + A2 * calmin * hV) / (A2 * hV + B2 * (1.0 - hV));
        ScalarValue c_beta = (A2 * (c - calmin * (1.0 - hV)) + B2 * cbtmin * (1.0 - hV)) /
                             (A2 * hV + B2 * (1.0 - hV));

        // Free energy for each phase and their first and second derivatives
        ScalarValue faV  = A0 + A2 * (c_alpha - calmin) * (c_alpha - calmin);
        ScalarValue fbV  = B0 + B2 * (c_beta - cbtmin) * (c_beta - cbtmin);
        ScalarValue fbcV = 2.0 * B2 * (c_beta - cbtmin);

        // Double-Well function (can be used to tune the interfacial energy)
        ScalarValue fbarriernV = 2.0 * n - 6.0 * n * n + 4.0 * n * n * n;

        // -------------------------------------------------
        // Nucleation expressions
        // -------------------------------------------------
        ScalarValue source_term(0.0);
        ScalarValue gamma(1.0);
        seed_nucleus(variable_list.get_q_point_location(), source_term, gamma, sim_timer);
        // -------------------------------------------------

        // Set the terms in the governing equations

        // For concentration
        ScalarValue eq_c = c;
        ScalarGrad  eqx_c =
          ScalarValue(-McV * delta_t) * (cx + (c_alpha - c_beta) * hnV * nx);

        // For order parameter (gamma is a variable order parameter mobility factor)
        ScalarValue eq_n = n - ScalarValue(delta_t * MnV) * gamma *
                                 ((fbV - faV) * hnV - (c_beta - c_alpha) * fbcV * hnV +
                                  W_barrier * fbarriernV);
        ScalarGrad eqx_n = ScalarValue(-delta_t * KnV * MnV) * gamma * nx;

        // --- Submitting the terms for the governing equations ---

        // Terms for the equation to evolve the concentration
        variable_list.set_value_term(0, eq_c);
        variable_list.set_gradient_term(0, eqx_c);

        // Terms for the equation to evolve the order parameter
        variable_list.set_value_term(1, eq_n + source_term);
        variable_list.set_gradient_term(1, eqx_n);
      }
    else if (solve_block_id == 1) // nucleation rate
      {
        ScalarValue c = variable_list.template get_value<Scalar, Current>(0);
        ScalarValue n = variable_list.template get_value<Scalar, Current>(1);
        // Supersaturation factor
        const ScalarValue ssf1 = c - calmin;
        ScalarValue       ssf(1.0);
        for (unsigned int d = 0; d < dim - 1; ++d)
          {
            ssf *= ssf1;
          }
        // Terms for the nucleation rate
        auto max = [](const ScalarValue &arr, number val)
        {
          ScalarValue result;
          for (unsigned int i = 0; i < arr.size(); ++i)
            {
              result[i] = std::max(arr[i], val);
            }
          return result;
        };
        // Calculate the nucleation rate
        double current_time = sim_timer.get_time();
        using std::exp;
        ScalarValue J =
          k1 * exp(-k2 / (max(ssf, number(1.0e-6)))) * exp(-tau / current_time);
        for (unsigned int i = 0; i < ScalarValue::size(); ++i)
          {
            if (n[i] >= 0.01)
              {
                J[i] = 0.0;
              }
          }
        variable_list.set_value_term(2, J);
      }
  }

  void
  seed_nucleus(const dealii::Point<dim, ScalarValue> &q_point_loc,
               ScalarValue                           &source_term,
               ScalarValue                           &gamma,
               const SimulationTimer                 &sim_timer) const
  {
    unsigned int current_increment = sim_timer.get_increment();
    double       current_time      = sim_timer.get_time();
    // Iterate through nuclei list
    for (const prisms::Nucleus<dim> &nucleus : get_pf_tools().nuclei_list)
      {
        // Calculate the distance function to the nucleus center
        const dealii::Point<dim, ScalarValue> loc_as_arr = [&]()
        {
          dealii::Point<dim, ScalarValue> result;
          const dealii::Point<dim>       &point = nucleus.location;
          for (unsigned int d = 0; d < dim; ++d)
            {
              result[d] = ScalarValue(point[d]);
            }
          return result;
        }();

        ScalarValue dist = distance<dim, ScalarValue>(
          q_point_loc,
          loc_as_arr,
          get_user_inputs().spatial_discretization.rectangular_mesh);
        // Seed a nucleus if it was added to the list of nuclei recently
        if (current_time < nucleus.seed_time + seeding_duration)
          {
            gamma *= 0.5 * (1.0 + std::tanh((dist - r_freeze) / interface_coeff));
          }
        if (nucleus.seed_increment == current_increment - 1)
          {
            source_term += 0.5 * (1.0 - std::tanh((dist - r_nuc) / interface_coeff));
          }
      }
  }
};

PRISMS_PF_END_NAMESPACE
