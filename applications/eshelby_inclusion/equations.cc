// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
CustomAttributeLoader::load_variable_attributes()
{
  set_variable_name(0, "u");
  set_variable_type(0, Vector);
  set_variable_equation_type(0, TimeIndependent);
  set_dependencies_gradient_term_rhs(0, "grad(u)");
  set_dependencies_gradient_term_lhs(0, "grad(change(u))");
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index solve_block) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index                                               solve_block,
  [[maybe_unused]] Types::Index current_index) const
{
  if (current_index == 0)
    {
      VectorGrad ux = variable_list.template get_symmetric_gradient<VectorGrad>(0);

      ScalarValue dist_from_inclusion = 0.0;
      ScalarValue inclusion_radius    = 10.0;
      for (unsigned int i = 0; i < dim; i++)
        {
          dist_from_inclusion += (q_point_loc[i] - 0.0) * (q_point_loc[i] - 0.0);
        }
      dist_from_inclusion = std::sqrt(dist_from_inclusion);

      VectorGrad transformation_strain;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              if (i == j)
                {
                  transformation_strain[i][j] =
                    0.01 * (0.5 + 0.5 *
                                    (1.0 - std::exp(-20.0 * (dist_from_inclusion -
                                                             inclusion_radius))) /
                                    (1.0 + std::exp(-20.0 * (dist_from_inclusion -
                                                             inclusion_radius))));
                }
              else
                {
                  transformation_strain[i][j] = 0.0;
                }
            }
        }
      VectorGrad stress;
      compute_stress<dim, ScalarValue>(compliance, ux - transformation_strain, stress);
      variable_list.set_gradient_term(0, -stress);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index                                               solve_block,
  [[maybe_unused]] Types::Index current_index) const
{
  if (current_index == 0)
    {
      VectorGrad change_ux =
        variable_list.template get_symmetric_gradient<VectorGrad>(0, Change);
      VectorGrad stress;
      compute_stress<dim, ScalarValue>(compliance, change_ux, stress);
      variable_list.set_gradient_term(0, stress, Change);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index solve_block) const
{}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
