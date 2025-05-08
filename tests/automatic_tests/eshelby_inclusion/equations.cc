// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

//JM added includes and prisms namespace
#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
customAttributeLoader::loadVariableAttributes()
{
  // Displacement Variable
  set_variable_name(0, "u");
  set_variable_type(0, VECTOR);
  set_variable_equation_type(0, TIME_INDEPENDENT);

  set_dependencies_gradient_term_RHS(0, "grad(u)");
  set_dependencies_gradient_term_LHS(0, "grad(change(u))");

  // Calculated Solution Variable (Probably unnecessary, but I'll fix syntax later)
  set_variable_name(1, "u_prisms");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(1, "u");
  set_is_postprocessed_field(1, true);

  // Difference Field Variable
  set_variable_name(2, "u_diff");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(2, "u");
  set_is_postprocessed_field(2, true);

  // Second Difference Field
  /*
  set_variable_name(4, "u_error");
  set_variable_type(4, SCALAR);
  set_variable_equation_type(4, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(4, "u");
  set_is_postprocessed_field(4, true);
  */
}

//JM updated template and customPDE
template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc) 
  const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) 
  const
{
  if (current_index == 0)
    {
      // u
      vectorGrad ux = variable_list.get_vector_symmetric_gradient(0);
      scalarValue dist_from_inclusion = 0.0;
      for (unsigned int i = 0; i < dim; i++)
        {
          dist_from_inclusion += (q_point_loc[i] - center[i]) *
                                 (q_point_loc[i] - center[i]);
        }
      dist_from_inclusion = std::sqrt(dist_from_inclusion);
      vectorGrad transformation_strain;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              if (i == j)
                {
                  for (unsigned int lane = 0; lane<dist_from_inclusion.size();lane++)
                    {
                      
                        transformation_strain[i][j][lane] = 
                          0.01 * (0.5 + 0.5 * (-1.0 * std::tanh(-20.0 * (dist_from_inclusion[lane] - incRadius))));
                    }
                }
              else
                {
                  transformation_strain[i][j] = 0.0;
                }
            }
        }
      vectorGrad stress;
      compute_stress<dim, scalarValue>(CIJ, ux - transformation_strain, stress);
      variable_list.set_vector_gradient_term(0, -stress);
    }
}

//JM helper function: a Kronecker Delta
double kDeltaEq(unsigned int i, unsigned int j) 
{
  if (i == j) 
    {
      return 1.0;
    } 
  else 
    {
      return 0.0;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_LHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>>  &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  if (current_index == 0)
    {
      vectorGrad change_ux = variable_list.get_vector_symmetric_gradient(0, CHANGE);
      vectorGrad stress;
      compute_stress<dim, scalarValue>(CIJ, change_ux, stress);
      variable_list.set_vector_gradient_term(0, stress, CHANGE);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_postprocess_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  // Adding in masking of inclusion
  scalarValue dist_from_inclusion = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      dist_from_inclusion += ((q_point_loc[i] - center[i]) * (q_point_loc[i] - center[i]));
    }

  for (unsigned int lane=0; lane<dist_from_inclusion.size(); lane++)
    {
      const double lane_distance = std::sqrt(dist_from_inclusion[lane]);
      if (incRadius >= lane_distance)
        {
          scalarValue u_prisms = 0.0;
          scalarValue diff = 0.0;

          variable_list.set_scalar_value_term(1, u_prisms);
          variable_list.set_scalar_value_term(2, diff);
        }
      else
        {
          vectorValue u_analytic_vector;
          vectorValue u_prisms_vector = variable_list.get_vector_value(0);
          scalarValue eshelbyConstant = (incRadius*incRadius*incRadius)/(6.0*(1-poisson));
          vectorGrad transformation_strain;
          for (unsigned int i = 0; i < dim; i++)
            {
              for (unsigned int j = 0; j < dim; j++)
                {
                  if (i == j)
                    {
                      for (unsigned int lane = 0; lane<dist_from_inclusion.size();lane++)
                        {
                          transformation_strain[i][j][lane] = 
                            0.01 * (0.5 + 0.5 * (-1.0 * std::tanh(-20.0 * (dist_from_inclusion[lane] - incRadius))));
                        }
                    }
                  else
                    {
                      transformation_strain[i][j] = 0.0;
                    }
                }
            }
          for (unsigned int i = 0; i < dim; i++)
            {
              double G = 0.0;
              for (unsigned int j = 0; j < dim; j ++)
                {
                  for (unsigned int k = 0; k < dim; k++)
                    {
                      if (j == k)
                        {
                          G += 0.01 * ((1.0 - 2.0 * poisson) *
                               (kDeltaEq(i,j) * ((q_point_loc[k][lane] - center[k]) / lane_distance) + 
                                kDeltaEq(i,k) * ((q_point_loc[j][lane] - center[j]) / lane_distance) - 
                               ((q_point_loc[i][lane] - center[i]) / lane_distance)) + 
                                3.0 * ((q_point_loc[i][lane] - center[i]) / lane_distance) * 
                               ((q_point_loc[j][lane] - center[j]) / lane_distance) * 
                               ((q_point_loc[k][lane] - center[k]) / lane_distance));
                        }
                      else
                        {
                          G += 0.0;
                        }
                    }
                }
              u_analytic_vector[i] = -1.0*eshelbyConstant*(1.0/(lane_distance*lane_distance))*G;
            }
          scalarValue u_analytic = u_analytic_vector.norm();
          scalarValue u_prisms = u_prisms_vector.norm();
          scalarValue diff_u = abs(u_analytic - u_prisms);

          variable_list.set_scalar_value_term(1, u_prisms);
          variable_list.set_scalar_value_term(2, diff_u);
        }
    }
  //dist_from_inclusion[] = std::sqrt(dist_from_inclusion);
}

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE