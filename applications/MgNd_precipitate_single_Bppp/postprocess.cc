// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

// =============================================================================================
// loadPostProcessorVariableAttributes: Set the attributes of the postprocessing
// variables
// =============================================================================================
// This function is analogous to 'load_variable_attributes' in 'equations.h', but
// for the postprocessing expressions. It sets the attributes for each
// postprocessing expression, including its name, whether it is a vector or
// scalar (only scalars are supported at present), its dependencies on other
// variables and their derivatives, and whether to calculate an integral of the
// postprocessed quantity over the entire domain. Note: this function is not a
// member of CustomPDE.

void
CustomAttributeLoader::loadPostProcessorVariableAttributes()
{
  // Variable 0
  set_variable_name(0, "f_tot");
  set_variable_type(0, FieldInfo::TensorRank::Scalar);

  set_dependencies_value_term_rhs(0, "c, grad(mu), n1, grad(n1), grad(u)");
  set_dependencies_gradient_term_rhs(0, "");

  set_output_integral(0, true);

  // Variable 1
  set_variable_name(1, "von_mises_stress");
  set_variable_type(1, FieldInfo::TensorRank::Scalar);

  set_dependencies_value_term_rhs(1, "c, n1, grad(u)");
  set_dependencies_gradient_term_rhs(1, "");

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
CustomPDE<dim, degree>::postProcessedFields(
  [[maybe_unused]] const VariableContainer<dim, degree, VectorizedArray<double>>
    &variable_list,
  [[maybe_unused]] VariableContainer<dim, degree, VectorizedArray<double>>
                                                            &pp_variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
  [[maybe_unused]] const VectorizedArray<double>             element_volume) const
{
  // --- Getting the values and derivatives of the model variables ---

  /// The concentration and its derivatives (names here should match those in
  /// the macros above)
  scalarvalueType c = variable_list.template get_value<ScalarValue>(0);

  // The first order parameter and its derivatives (names here should match
  // those in the macros above)
  scalarvalueType n1  = variable_list.template get_value<ScalarValue>(2);
  scalargradType  n1x = variable_list.template get_gradient<ScalarGrad>(2);

  // The derivative of the displacement vector (names here should match those in
  // the macros above)
  vectorgradType ux = variable_list.template get_gradient<VectorGrad>(3);

  // --- Setting the expressions for the terms in the postprocessing expressions
  // ---

  scalarvalueType h1V  = (3.0 * n1 * n1 - 2.0 * n1 * n1 * n1);
  scalarvalueType hn1V = (6.0 * n1 - 6.0 * n1 * n1);

  // This double-well function can be used to tune the interfacial energy
  scalarvalueType fbarrierV = (n1 * n1 - 2.0 * n1 * n1 * n1 + n1 * n1 * n1 * n1);

  // Calculate c_alpha and c_beta from c
  scalarvalueType c_alpha =
    ((B2 * c + 0.5 * (B1 - A1) * h1V) / (A2 * h1V + B2 * (1.0 - h1V)));
  scalarvalueType c_beta =
    ((A2 * c + 0.5 * (A1 - B1) * (1.0 - h1V)) / (A2 * h1V + B2 * (1.0 - h1V)));

  scalarvalueType faV   = (A2 * c_alpha * c_alpha + A1 * c_alpha + A0);
  scalarvalueType facV  = (2.0 * A2 * c_alpha + A1);
  scalarvalueType faccV = (constV(2.0) * A2);
  scalarvalueType fbV   = (B2 * c_beta * c_beta + B1 * c_beta + B0);
  scalarvalueType fbcV  = (2.0 * B2 * c_beta + B1);
  scalarvalueType fbccV = (constV(2.0) * B2);

  // Start calculating components of the energy density
  scalarvalueType total_energy_density = constV(0.0);

  scalarvalueType f_chem =
    (constV(1.0) - (h1V)) * faV + (h1V) *fbV + constV(W) * fbarrierV;

  scalarvalueType f_grad = constV(0.0);

  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += constV(0.5 * Kn1[i][j]) * n1x[i] * n1x[j];
        }
    }

  // Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
  scalarvalueType cbnV, cbcV, cbcnV, cacV;

  cbcV  = faccV / ((constV(1.0) - h1V) * fbccV + h1V * faccV);
  cacV  = fbccV / ((constV(1.0) - h1V) * fbccV + h1V * faccV);
  cbnV  = hn1V * (c_alpha - c_beta) * cbcV;
  cbcnV = (faccV * (fbccV - faccV) * hn1V) /
          (((1.0 - h1V) * fbccV + h1V * faccV) *
           ((1.0 - h1V) * fbccV +
            h1V * faccV)); // Note: this is only true if faV and fbV are quadratic

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  Tensor<2, dim, VectorizedArray<double>> sfts1, sfts1c, sfts1cc, sfts1n, sfts1cn;

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c_beta + b_p
          sfts1[i][j]   = constV(sfts_linear1[i][j]) * c_beta + constV(sfts_const1[i][j]);
          sfts1c[i][j]  = constV(sfts_linear1[i][j]) * cbcV;
          sfts1cc[i][j] = constV(0.0);
          sfts1n[i][j]  = constV(sfts_linear1[i][j]) * cbnV;
          sfts1cn[i][j] = constV(sfts_linear1[i][j]) * cbcnV;
        }
    }

  // compute E2=(E-E0)
  VectorizedArray<double> E2[dim][dim], S[dim][dim];

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E2[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]) - (sfts1[i][j] * h1V);
        }
    }

  // compute stress
  // S=C*(E-E0)
  VectorizedArray<double> CIJ_combined[2 * dim - 1 + dim / 3][2 * dim - 1 + dim / 3];

  if (n_dependent_stiffness == true)
    {
      for (unsigned int i = 0; i < 2 * dim - 1 + dim / 3; i++)
        {
          for (unsigned int j = 0; j < 2 * dim - 1 + dim / 3; j++)
            {
              CIJ_combined[i][j] =
                CIJ_Mg[i][j] * (constV(1.0) - h1V) + CIJ_Beta[i][j] * h1V;
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

  total_energy_density = f_chem + f_grad + f_el;

  // Calculate the chemical potential for the concentration
  scalarvalueType mu_c = constV(0.0);
  mu_c += facV * cacV * (constV(1.0) - h1V) + fbcV * cbcV * h1V;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          mu_c -= S[i][j] * (sfts1c[i][j] * h1V);
        }
    }

  scalarvalueType mu_c_el = constV(0.0);
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          mu_c_el -= S[i][j] * (sfts1c[i][j] * h1V);
        }
    }

  // The Von Mises Stress
  VectorizedArray<double> vm_stress;
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

  // Compute one of the stress terms in the order parameter chemical potential,
  // nDependentMisfitACp = -C*(E-E0)*(E0_n)
  VectorizedArray<double> nDependentMisfitAC1    = constV(0.0);
  VectorizedArray<double> nDependentMisfitAC1_t1 = constV(0.0);
  VectorizedArray<double> nDependentMisfitAC1_t2 = constV(0.0);

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          nDependentMisfitAC1 += -S[i][j] * (sfts1n[i][j] * h1V + sfts1[i][j] * hn1V);
          nDependentMisfitAC1_t1 += -S[i][j] * (sfts1n[i][j] * h1V);
          nDependentMisfitAC1_t2 += -S[i][j] * (sfts1[i][j] * hn1V);
        }
    }

  // Compute the other stress term in the order parameter chemical potential,
  // heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
  VectorizedArray<double> heterMechAC1 = constV(0.0);
  VectorizedArray<double> S2[dim][dim];

  if (n_dependent_stiffness == true)
    {
      computeStress<dim>(CIJ_Beta - CIJ_Mg, E2, S2);

      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              heterMechAC1 += S2[i][j] * E2[i][j];
            }
        }
      heterMechAC1 = 0.5 * hn1V * heterMechAC1;
    }

  // compute K*nx
  scalargradType Knx1;
  for (unsigned int a = 0; a < dim; a++)
    {
      Knx1[a] = 0.0;
      for (unsigned int b = 0; b < dim; b++)
        {
          Knx1[a] += constV(Kn1[a][b]) * n1x[b];
        }
    }

  // --- Submitting the terms for the postprocessing expressions ---

  pp_variable_list.set_scalar_value_term_rhs(0, total_energy_density);
  pp_variable_list.set_scalar_value_term_rhs(1, vm_stress);
}