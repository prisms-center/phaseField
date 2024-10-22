// =============================================================================================
// loadPostProcessorVariableAttributes: Set the attributes of the postprocessing
// variables
// =============================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but
// for the postprocessing expressions. It sets the attributes for each
// postprocessing expression, including its name, whether it is a vector or
// scalar (only scalars are supported at present), its dependencies on other
// variables and their derivatives, and whether to calculate an integral of the
// postprocessed quantity over the entire domain. Note: this function is not a
// member of customPDE.

void
variableAttributeLoader::loadPostProcessorVariableAttributes()
{
  // Variable 0
  set_variable_name(0, "f_tot");
  set_variable_type(0, SCALAR);

  set_dependencies_value_term_RHS(
    0,
    "c, grad(c), n1, grad(n1), n2, grad(n2), n3, grad(n3), grad(u)");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(0, true);
}

// =================================================================================
// Define the expressions for the post-processed fields
// =================================================================================

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
  // --- Getting the values and derivatives of the model variables ---

  // The concentration and its derivatives
  scalarvalueType c = variable_list.get_scalar_value(0);

  // The first order parameter and its derivatives
  scalarvalueType n1  = variable_list.get_scalar_value(1);
  scalargradType  n1x = variable_list.get_scalar_gradient(1);

  // The second order parameter and its derivatives
  scalarvalueType n2  = variable_list.get_scalar_value(2);
  scalargradType  n2x = variable_list.get_scalar_gradient(2);

  // The third order parameter and its derivatives
  scalarvalueType n3  = variable_list.get_scalar_value(3);
  scalargradType  n3x = variable_list.get_scalar_gradient(3);

  // The derivative of the displacement vector
  vectorgradType ux = variable_list.get_vector_gradient(4);

  // --- Setting the expressions for the terms in the postprocessing expressions
  // ---

  scalarvalueType f_tot = constV(0.0);

  // Free energy expressions and interpolation functions
  scalarvalueType faV = pfunct_faV.val(c);
  scalarvalueType fbV = pfunct_fbV.val(c);
  scalarvalueType h1V =
    (10.0 * n1 * n1 * n1 - 15.0 * n1 * n1 * n1 * n1 + 6.0 * n1 * n1 * n1 * n1 * n1);
  scalarvalueType h2V =
    (10.0 * n2 * n2 * n2 - 15.0 * n2 * n2 * n2 * n2 + 6.0 * n2 * n2 * n2 * n2 * n2);
  scalarvalueType h3V =
    (10.0 * n3 * n3 * n3 - 15.0 * n3 * n3 * n3 * n3 + 6.0 * n3 * n3 * n3 * n3 * n3);

  scalarvalueType f_chem =
    (constV(1.0) - (h1V + h2V + h3V)) * faV + (h1V + h2V + h3V) * fbV;

  scalarvalueType f_grad = constV(0.0);

  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += constV(0.5 * Kn1[i][j]) * n1x[i] * n1x[j];
        }
    }

  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += constV(0.5 * Kn2[i][j]) * n2x[i] * n2x[j];
        }
    }

  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += constV(0.5 * Kn3[i][j]) * n3x[i] * n3x[j];
        }
    }

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  Tensor<2, dim, VectorizedArray<double>> sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc,
    sfts3, sfts3c, sfts3cc;

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts1[i][j]   = constV(sfts_linear1[i][j]) * c + constV(sfts_const1[i][j]);
          sfts1c[i][j]  = constV(sfts_linear1[i][j]);
          sfts1cc[i][j] = constV(0.0);

          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts2[i][j]   = constV(sfts_linear2[i][j]) * c + constV(sfts_const2[i][j]);
          sfts2c[i][j]  = constV(sfts_linear1[i][j]);
          sfts2cc[i][j] = constV(0.0);

          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts3[i][j]   = constV(sfts_linear3[i][j]) * c + constV(sfts_const3[i][j]);
          sfts3c[i][j]  = constV(sfts_linear3[i][j]);
          sfts3cc[i][j] = constV(0.0);
        }
    }

  // compute E2=(E-E0)
  VectorizedArray<double> E2[dim][dim], S[dim][dim];

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          // E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V +
          // sfts2[i][j]*h2V + sfts3[i][j]*h3V);
          E2[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]) -
                     (sfts1[i][j] * h1V + sfts2[i][j] * h2V + sfts3[i][j] * h3V);
        }
    }

  // compute stress
  // S=C*(E-E0)
  VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

  if (n_dependent_stiffness == true)
    {
      VectorizedArray<double> sum_hV;
      sum_hV = h1V + h2V + h3V;
      for (unsigned int i = 0; i < 2 * dim - 1 + dim / 3; i++)
        {
          for (unsigned int j = 0; j < 2 * dim - 1 + dim / 3; j++)
            {
              CIJ_combined[i][j] =
                CIJ_Mg[i][j] * (constV(1.0) - sum_hV) + CIJ_Beta[i][j] * sum_hV;
            }
        }
      computeStress<dim>(CIJ_combined, E2, S);
    }
  else
    {
      computeStress<dim>(CIJ_Mg, E2, S);
    }

  scalarvalueType f_el = constV(0.0);

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_el += constV(0.5) * S[i][j] * E2[i][j];
        }
    }

  f_tot = f_chem + f_grad + f_el;

  // --- Submitting the terms for the postprocessing expressions ---

  pp_variable_list.set_scalar_value_term_RHS(0, f_tot);
}
