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

  set_dependencies_value_term_RHS(0,
                                  "c, n1, grad(n1), n2, grad(n2), n2, grad(n3), grad(u)");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(0, true);

  // Variable 1
  set_variable_name(1, "f_el");
  set_variable_type(1, SCALAR);

  set_dependencies_value_term_RHS(0, "grad(u)");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(1, false);

  // Variable 2
  set_variable_name(2, "von_mises_stress");
  set_variable_type(2, SCALAR);

  set_dependencies_value_term_RHS(0, "grad(u)");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(2, false);
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

  // The interpolation functions and their derivatives
  scalarvalueType h1V  = (n1 * n1 * n1 * (10.0 - 15.0 * n1 + 6.0 * n1 * n1));
  scalarvalueType hn1V = (30.0 * (n1 - 1.0) * (n1 - 1.0) * n1 * n1);
  scalarvalueType h2V  = (n2 * n2 * n2 * (10.0 - 15.0 * n2 + 6.0 * n2 * n2));
  scalarvalueType hn2V = (30.0 * (n2 - 1.0) * (n2 - 1.0) * n2 * n2);
  scalarvalueType h3V  = (n3 * n3 * n3 * (10.0 - 15.0 * n3 + 6.0 * n3 * n3));
  scalarvalueType hn3V = (30.0 * (n3 - 1.0) * (n3 - 1.0) * n3 * n3);

  scalarvalueType sum_hpV = h1V + h2V + h3V;
  scalarvalueType c_alpha =
    ((B2 * c + 0.5 * (B1 - A1) * sum_hpV)) / (A2 * (sum_hpV) + B2 * (1.0 - sum_hpV));
  scalarvalueType c_beta = ((A2 * c + 0.5 * (A1 - B1) * (1.0 - sum_hpV)) /
                            (A2 * (sum_hpV) + B2 * (1.0 - sum_hpV)));

  // The free energy expressions and their derivatives
  scalarvalueType faV   = (A2 * c_alpha * c_alpha + A1 * c_alpha + A0);
  scalarvalueType facV  = (2.0 * A2 * c_alpha + A1);
  scalarvalueType faccV = (constV(2.0) * A2);
  scalarvalueType fbV   = (B2 * c_beta * c_beta + B1 * c_beta + B0);
  scalarvalueType fbcV  = (2.0 * B2 * c_beta + B1);
  scalarvalueType fbccV = (constV(2.0) * B2);

  // This double-well function can be used to tune the interfacial energy and
  // its derivatives
  scalarvalueType fbarrierV =
    (n1 * n1 - 2.0 * n1 * n1 * n1 + n1 * n1 * n1 * n1 + n2 * n2 - 2.0 * n2 * n2 * n2 +
     n2 * n2 * n2 * n2 + n3 * n3 - 2.0 * n3 * n3 * n3 + n3 * n3 * n3 * n3 +
     5.0 * (n1 * n1 * n2 * n2 + n1 * n1 * n3 * n3 + n2 * n2 * n3 * n3) +
     5.0 * n1 * n1 * n2 * n2 * n3 * n3);
  scalarvalueType fbarriern1V =
    (2.0 * n1 - 6.0 * n1 * n1 + 4.0 * n1 * n1 * n1 + 10.0 * n1 * (n2 * n2 + n3 * n3) +
     10.0 * n1 * n2 * n2 * n3 * n3);
  scalarvalueType fbarriern2V =
    (2.0 * n2 - 6.0 * n2 * n2 + 4.0 * n2 * n2 * n2 + 10.0 * n2 * (n1 * n1 + n3 * n3) +
     10.0 * n2 * n1 * n1 * n3 * n3);
  scalarvalueType fbarriern3V =
    (2.0 * n3 - 6.0 * n3 * n3 + 4.0 * n3 * n3 * n3 + 10.0 * n3 * (n2 * n2 + n1 * n1) +
     10.0 * n3 * n2 * n2 * n1 * n1);

  // Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
  // Note: this section can be optimized to reduce recalculations
  scalarvalueType cbn1V, cbn2V, cbn3V, cbcV, cbcn1V, cbcn2V, cbcn3V, cacV;

  cbcV = faccV / ((constV(1.0) - sum_hpV) * fbccV + (sum_hpV) *faccV);
  cacV = fbccV / ((constV(1.0) - sum_hpV) * fbccV + (sum_hpV) *faccV);

  cbn1V = hn1V * (c_alpha - c_beta) * cbcV;
  cbn2V = hn2V * (c_alpha - c_beta) * cbcV;
  cbn3V = hn3V * (c_alpha - c_beta) * cbcV;

  cbcn1V = (faccV * (fbccV - faccV) * hn1V) /
           (((1.0 - sum_hpV) * fbccV + sum_hpV * faccV) *
            ((1.0 - sum_hpV) * fbccV +
             sum_hpV * faccV)); // Note: this is only true if faV and fbV are quadratic
  cbcn2V = (faccV * (fbccV - faccV) * hn2V) /
           (((1.0 - sum_hpV) * fbccV + sum_hpV * faccV) *
            ((1.0 - sum_hpV) * fbccV +
             sum_hpV * faccV)); // Note: this is only true if faV and fbV are quadratic
  cbcn3V = (faccV * (fbccV - faccV) * hn3V) /
           (((1.0 - sum_hpV) * fbccV + sum_hpV * faccV) *
            ((1.0 - sum_hpV) * fbccV +
             sum_hpV * faccV)); // Note: this is only true if faV and fbV are quadratic

  scalarvalueType f_chem =
    (constV(1.0) - (h1V + h2V + h3V)) * faV + (h1V + h2V + h3V) * fbV + W * (fbarrierV);

  scalarvalueType f_grad = constV(0.0);

  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += constV(0.5 * Kn1[i][j]) * n1x[i] * n1x[j] +
                    constV(0.5 * Kn2[i][j]) * n2x[i] * n2x[j] +
                    constV(0.5 * Kn3[i][j]) * n3x[i] * n3x[j];
        }
    }

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>> sfts1, sfts1c, sfts1cc, sfts2,
    sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c_beta + b_p
          sfts1[i][j] = constV(sfts_const1[i][j]);
          sfts2[i][j] = constV(sfts_const2[i][j]);
          sfts3[i][j] = constV(sfts_const3[i][j]);
        }
    }

  // compute E2=(E-E0)
  scalarvalueType E2[dim][dim], S[dim][dim];

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E2[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]) -
                     (sfts1[i][j] * h1V + sfts2[i][j] * h2V + sfts3[i][j] * h3V);
        }
    }

  // compute stress
  // S=C*(E-E0)
  scalarvalueType CIJ_combined[2 * dim - 1 + dim / 3][2 * dim - 1 + dim / 3];

  if (n_dependent_stiffness == true)
    {
      for (unsigned int i = 0; i < 2 * dim - 1 + dim / 3; i++)
        {
          for (unsigned int j = 0; j < 2 * dim - 1 + dim / 3; j++)
            {
              CIJ_combined[i][j] = CIJ_Mg[i][j] * (constV(1.0) - sum_hpV) +
                                   CIJ_Beta1[i][j] * (h1V) + CIJ_Beta2[i][j] * (h2V) +
                                   CIJ_Beta3[i][j] * (h3V);
            }
        }
      computeStress<dim>(CIJ_combined, E2, S);
    }
  else
    {
      computeStress<dim>(CIJ_Mg, E2, S);
    }

  dealii::VectorizedArray<double> f_el = constV(0.0);

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_el += constV(0.5) * S[i][j] * E2[i][j];
        }
    }

  scalarvalueType total_energy_density = f_chem + f_grad + f_el;

  dealii::VectorizedArray<double> vm_stress;
  if (dim == 3)
    {
      vm_stress = (S[0][0] - S[1][1]) * (S[0][0] - S[1][1]) +
                  (S[1][1] - S[2][2]) * (S[1][1] - S[2][2]) +
                  (S[2][2] - S[0][0]) * (S[2][2] - S[0][0]);
      vm_stress +=
        constV(6.0) * (S[0][1] * S[0][1] + S[1][2] * S[1][2] + S[2][0] * S[2][0]);
      vm_stress *= constV(0.5);
      vm_stress = std::sqrt(vm_stress);
    }
  else
    {
      vm_stress = S[0][0] * S[0][0] - S[0][0] * S[1][1] + S[1][1] * S[1][1] +
                  constV(3.0) * S[0][1] * S[0][1];
      vm_stress = std::sqrt(vm_stress);
    }

  // --- Submitting the terms for the postprocessing expressions ---
  pp_variable_list.set_scalar_value_term_RHS(0, total_energy_density);
  pp_variable_list.set_scalar_value_term_RHS(1, f_el);
  pp_variable_list.set_scalar_value_term_RHS(2, vm_stress);
}
