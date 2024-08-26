// =================================================================================
// Set the attributes of the postprocessing variables
// =================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but
// for the postprocessing expressions. It sets the attributes for each
// postprocessing expression, including its name, whether it is a vector or
// scalar (only scalars are supported at present), its dependencies on other
// variables and their derivatives, and whether to calculate an integral of the
// postprocessed quantity over the entire domain.

void
variableAttributeLoader::loadPostProcessorVariableAttributes()
{
  // Variable 0
  set_variable_name(0, "feature_ids");
  set_variable_type(0, SCALAR);

  set_dependencies_value_term_RHS(0, "n0, n1, n2, n3, n4, n5");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(0, false);

  // Variable 1
  set_variable_name(1, "op_ids");
  set_variable_type(1, SCALAR);

  set_dependencies_value_term_RHS(1, "n0, n1, n2, n3, n4, n5");
  set_dependencies_gradient_term_RHS(1, "");

  set_output_integral(1, false);
}

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================
// This function is analogous to 'explicitEquationRHS' and
// 'nonExplicitEquationRHS' in equations.h. It takes in "variable_list" and
// "q_point_loc" as inputs and outputs two terms in the expression for the
// postprocessing variable -- one proportional to the test function and one
// proportional to the gradient of the test function. The index for each
// variable in this list corresponds to the index given at the top of this file
// (for submitting the terms) and the index in 'equations.h' for assigning the
// values/derivatives of the primary variables.

template <int dim, int degree>
void
customPDE<dim, degree>::postProcessedFields(
  const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  variableContainer<dim, degree, dealii::VectorizedArray<double>>       &pp_variable_list,
  const dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const
{
  // --- Getting the values and derivatives of the model variables ---

  scalarvalueType ni;

  scalarvalueType max_val = constV(-1.0);
  scalarvalueType max_op  = constV(100.0);
  for (unsigned int i = 0; i < userInputs.number_of_variables; i++)
    {
      ni = variable_list.get_scalar_value(i);

      for (unsigned int v = 0; v < ni.size(); v++)
        {
          if (ni[v] > max_val[v])
            {
              max_val[v] = ni[v];
              max_op[v]  = i;
            }
        }
    }

  scalarvalueType feature_ids = constV(-1.0);
  for (unsigned int v = 0; v < ni.size(); v++)
    {
      for (unsigned int g = 0; g < this->simplified_grain_representations.size(); g++)
        {
          unsigned int max_op_nonvec = (unsigned int) std::abs(max_op[v]);

          if (this->simplified_grain_representations[g].getOrderParameterId() ==
              max_op_nonvec)
            {
              dealii::Point<dim> q_point_loc_nonvec;
              for (unsigned int d = 0; d < dim; d++)
                {
                  q_point_loc_nonvec(d) = q_point_loc(d)[v];
                }

              double dist = 0.0;
              for (unsigned int d = 0; d < dim; d++)
                {
                  dist += (q_point_loc_nonvec(d) -
                           this->simplified_grain_representations[g].getCenter()(d)) *
                          (q_point_loc_nonvec(d) -
                           this->simplified_grain_representations[g].getCenter()(d));
                }
              dist = std::sqrt(dist);

              if (dist < (this->simplified_grain_representations[g].getRadius() +
                          userInputs.buffer_between_grains / 2.0))
                {
                  feature_ids[v] =
                    (double) (this->simplified_grain_representations[g].getGrainId());
                }
            }
        }
    }

  scalarvalueType sum_n = constV(0.0);
  for (unsigned int i = 0; i < userInputs.number_of_variables; i++)
    {
      ni = variable_list.get_scalar_value(i);
      sum_n += ni;
    }
  for (unsigned int v = 0; v < ni.size(); v++)
    {
      if (sum_n[v] < 0.01)
        {
          max_op[v]      = -1.0;
          feature_ids[v] = -1.0;
        }
    }

  // --- Submitting the terms for the postprocessing expressions ---

  pp_variable_list.set_scalar_value_term_RHS(0, feature_ids);
  pp_variable_list.set_scalar_value_term_RHS(1, max_op);
}
