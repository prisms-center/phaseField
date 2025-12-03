// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <deal.II/base/vectorization.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/config.h>
#include <prismspf/nucleation/nucleus.h>

#include "prismspf/utilities/periodic_distance.h"

PRISMS_PF_BEGIN_NAMESPACE

void
CustomAttributeLoader::load_variable_attributes()
{
  // Variable 0
  set_variable_name(0, "c");
  set_variable_type(0, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(0, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(0, "c");
  set_dependencies_gradient_term_rhs(0, "c, grad(c), n, grad(n)");

  // Variable 1
  set_variable_name(1, "n");
  set_variable_type(1, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(1, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(1, "c, n");
  set_dependencies_gradient_term_rhs(1, "grad(n)");

  // Variable 2
  set_variable_name(2, "nucleation_rate");
  set_variable_type(2, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(2, ExplicitTimeDependent);

  insert_dependencies_value_term_rhs(2, std::set<std::string> {"c", "n"});
  set_is_nucleation_rate(2, true, std::set<std::string> {"n"});
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  // --- Getting the values and derivatives of the model variables ---

  // The concentration and its derivatives
  ScalarValue c  = variable_list.template get_value<ScalarValue>(0);
  ScalarGrad  cx = variable_list.template get_gradient<ScalarGrad>(0);

  // The order parameter and its derivatives
  ScalarValue n  = variable_list.template get_value<ScalarValue>(1);
  ScalarGrad  nx = variable_list.template get_gradient<ScalarGrad>(1);

  // --- Setting the expressions for the terms in the governing equations ---

  // Interpolation function and its derivative
  ScalarValue hV  = 3.0 * n * n - 2.0 * n * n * n;
  ScalarValue hnV = 6.0 * n - 6.0 * n * n;

  double delta_t = get_timestep();

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
  ScalarValue source_term = ScalarValue(0.0);
  ScalarValue gamma       = ScalarValue(1.0);
  seed_nucleus(q_point_loc, source_term, gamma);
  // -------------------------------------------------

  // Set the terms in the governing equations

  // For concentration
  ScalarValue eq_c  = c;
  ScalarGrad  eqx_c = ScalarValue(-McV * delta_t) * (cx + (c_alpha - c_beta) * hnV * nx);

  // For order parameter (gamma is a variable order parameter mobility factor)
  ScalarValue eq_n =
    n - ScalarValue(delta_t * MnV) * gamma *
          ((fbV - faV) * hnV - (c_beta - c_alpha) * fbcV * hnV + W_barrier * fbarriernV);
  ScalarGrad eqx_n = ScalarValue(-delta_t * KnV * MnV) * gamma * nx;

  // Supersaturation factor
  const ScalarValue ssf1 = c - calmin;
  ScalarValue       ssf(1.0);
  for (unsigned int d = 0; d < dim - 1; ++d)
    {
      ssf *= ssf1;
    }

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
  double current_time = get_user_inputs().get_temporal_discretization().get_time();
  using std::exp;
  ScalarValue J = k1 * exp(-k2 / (max(ssf, number(1.0e-6)))) * exp(-tau / current_time);
  for (unsigned int i = 0; i < ScalarValue::size(); ++i)
    {
      if (n[i] >= 0.01)
        {
          J[i] = 0.0;
        }
    }

  // --- Submitting the terms for the governing equations ---

  // Terms for the equation to evolve the concentration
  variable_list.set_value_term(0, eq_c);
  variable_list.set_gradient_term(0, eqx_c);

  // Terms for the equation to evolve the order parameter
  variable_list.set_value_term(1, eq_n + source_term);
  variable_list.set_gradient_term(1, eqx_n);

  // Terms for the nucleation rate
  variable_list.set_value_term(2, J);
}

// =================================================================================
// seedNucleus: a function particular to this app
// =================================================================================

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::seed_nucleus(
  const dealii::Point<dim, ScalarValue> &q_point_loc,
  ScalarValue                           &source_term,
  ScalarValue                           &gamma) const
{
  unsigned int current_increment =
    get_user_inputs().get_temporal_discretization().get_increment();
  double current_time = get_user_inputs().get_temporal_discretization().get_time();
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

      ScalarValue dist =
        prisms::distance<dim, ScalarValue>(q_point_loc, loc_as_arr, get_user_inputs());
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

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           current_index) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           current_index) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE