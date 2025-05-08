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
customAttributeLoader::loadPostProcessorVariableAttributes()
{
  // TODO (Phil): enable submission of the following several values,
  //               via user input parameters

  // Number of order parameters
  unsigned int number_ops_ini = 23;

  // Number of empty OPs (for nucleation)
  unsigned int number_ops_nuc = 60;

  // Number of dislocation density fields
  unsigned int N_rho = 83;

  // Total number of order parameters
  unsigned int number_ops_total = number_ops_ini + number_ops_nuc;

  std::string all_op_names = "";
  for (unsigned int i = 0; i < number_ops_total; i++)
    {
      all_op_names.append("n" + std::to_string(i));
      if (i < number_ops_total - 1)
        {
          all_op_names.append(", ");
        }
    }

  // Variable 0
  set_variable_name(0, "feature_ids");
  set_variable_type(0, SCALAR);

  set_dependencies_value_term_RHS(0, all_op_names);
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(0, false);

  // Variable 1
  set_variable_name(1, "op_ids");
  set_variable_type(1, SCALAR);

  set_dependencies_value_term_RHS(1, all_op_names);
  set_dependencies_gradient_term_RHS(1, "");

  set_output_integral(1, false);

  // Variable 2
  set_variable_name(2, "sum_sq_op_ids");
  set_variable_type(2, SCALAR);

  set_dependencies_value_term_RHS(2, all_op_names);
  set_dependencies_gradient_term_RHS(2, "");

  set_output_integral(2, false);

  // Variable 3
  set_variable_name(3, "sum_rec_ids");
  set_variable_type(3, SCALAR);

  set_dependencies_value_term_RHS(3, all_op_names);
  set_dependencies_gradient_term_RHS(3, "");

  set_output_integral(3, false);
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
  [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
    &variable_list,
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                            &pp_variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
  [[maybe_unused]] const VectorizedArray<double>             element_volume) const
{
  // TODO (Phil): enable submission of the following several values,
  //               via user input parameters

  // Number of order parameters
  unsigned int number_ops_ini = 23;

  // Number of empty OPs (for nucleation)
  unsigned int number_ops_nuc = 60;

  // Number of dislocation density fields
  unsigned int N_rho = 83;

  // Total number of order parameters
  unsigned int number_ops_total = number_ops_ini + number_ops_nuc;

  // --- Getting the values and derivatives of the model variables ---

  scalarvalueType ni;

  scalarvalueType max_val = constV(-100.0);
  scalarvalueType max_op  = constV(100.0);
  for (unsigned int i = 0; i < number_ops_total; i++)
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
              Point<dim> q_point_loc_nonvec;
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
  for (unsigned int i = 0; i < number_ops_total; i++)
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

  scalarvalueType sum_n_sq = constV(0.0);
  for (unsigned int i = 0; i < number_ops_total; i++)
    {
      ni = variable_list.get_scalar_value(i);
      sum_n_sq += ni * ni;
    }

  scalarvalueType sum_n_nuc = constV(0.0);
  for (unsigned int i = number_ops_ini; i < number_ops_total; i++)
    {
      ni = variable_list.get_scalar_value(i);
      sum_n_nuc += ni;
    }

  // --- Submitting the terms for the postprocessing expressions ---

  pp_variable_list.set_scalar_value_term_RHS(0, feature_ids);
  pp_variable_list.set_scalar_value_term_RHS(1, max_op);
  pp_variable_list.set_scalar_value_term_RHS(2, sum_n_sq);
  pp_variable_list.set_scalar_value_term_RHS(3, sum_n_nuc);
}
