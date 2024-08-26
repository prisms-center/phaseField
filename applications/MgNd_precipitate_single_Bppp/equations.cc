// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for
// each function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void
variableAttributeLoader::loadVariableAttributes()
{
  // Variable 0
  set_variable_name(0, "c");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "c");
  set_dependencies_gradient_term_RHS(0, "n1, grad(mu)");

  // Variable 0
  set_variable_name(1, "mu");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, AUXILIARY);

  set_dependencies_value_term_RHS(1, "c, n1, grad(u)");
  set_dependencies_gradient_term_RHS(1, "");

  // Variable 1
  set_variable_name(2, "n1");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(2, "c, n1, grad(u)");
  set_dependencies_gradient_term_RHS(2, "grad(n1)");

  // Variable 2
  set_variable_name(3, "u");
  set_variable_type(3, VECTOR);
  set_variable_equation_type(3, TIME_INDEPENDENT);

  set_dependencies_value_term_RHS(3, "");
  set_dependencies_gradient_term_RHS(3, "n1, grad(u)");
  set_dependencies_value_term_LHS(3, "");
  set_dependencies_gradient_term_LHS(3, "n1, grad(change(u))");
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time
// dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a
// list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one
// proportional to the test function and one proportional to the gradient of the
// test function. The index for each variable in this list corresponds to the
// index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const
{
  // --- Getting the values and derivatives of the model variables ---

  // The concentration and its derivatives (names here should match those in the
  // macros above)
  scalarvalueType c = variable_list.get_scalar_value(0);

  scalargradType mux = variable_list.get_scalar_gradient(1);

  // The first order parameter and its derivatives (names here should match
  // those in the macros above)
  scalarvalueType n1  = variable_list.get_scalar_value(2);
  scalargradType  n1x = variable_list.get_scalar_gradient(2);

  // The derivative of the displacement vector (names here should match those in
  // the macros above)
  vectorgradType ux = variable_list.get_vector_gradient(3);

  // --- Setting the expressions for the terms in the governing equations ---

  scalarvalueType h1V  = (3.0 * n1 * n1 - 2.0 * n1 * n1 * n1);
  scalarvalueType hn1V = (6.0 * n1 - 6.0 * n1 * n1);

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

  // This double-well function can be used to tune the interfacial energy
  scalarvalueType fbarrierV  = (n1 * n1 - 2.0 * n1 * n1 * n1 + n1 * n1 * n1 * n1);
  scalarvalueType fbarriernV = (2.0 * n1 - 6.0 * n1 * n1 + 4.0 * n1 * n1 * n1);

  // Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
  scalarvalueType cacV, canV, cbnV, cbcV, cbcnV;

  cacV = fbccV / ((constV(1.0) - h1V) * fbccV + h1V * faccV);
  canV = hn1V * (c_alpha - c_beta) * cacV;

  cbcV  = faccV / ((constV(1.0) - h1V) * fbccV + h1V * faccV);
  cbnV  = hn1V * (c_alpha - c_beta) * cbcV;
  cbcnV = (faccV * (fbccV - faccV) * hn1V) /
          (((1.0 - h1V) * fbccV + h1V * faccV) *
           ((1.0 - h1V) * fbccV +
            h1V * faccV)); // Note: this is only true if faV and fbV are quadratic

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>> sfts1, sfts1c, sfts1cc, sfts1n,
    sfts1cn;

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
  dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E2[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]) - (sfts1[i][j] * h1V);
        }
    }

  // compute stress
  // S=C*(E-E0)
  //  Compute stress tensor (which is equal to the residual, Rux)
  dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

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

  // Compute one of the stress terms in the order parameter chemical potential,
  // nDependentMisfitACp = -C*(E-E0)*(E0_n)
  dealii::VectorizedArray<double> nDependentMisfitAC1 = constV(0.0);

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          nDependentMisfitAC1 += -S[i][j] * (sfts1n[i][j] * h1V + sfts1[i][j] * hn1V);
        }
    }

  // Compute the other stress term in the order parameter chemical potential,
  // heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
  dealii::VectorizedArray<double> heterMechAC1 = constV(0.0);
  dealii::VectorizedArray<double> S2[dim][dim];

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

  // The terms in the governing equations
  scalarvalueType eq_c = (c);
  scalargradType  eqx_c =
    (constV(-userInputs.dtValue * McV) * (h1V * faccV + (constV(1.0) - h1V) * fbccV) /
     (faccV * fbccV) * mux);

  scalarvalueType eq_n1  = (n1 - constV(userInputs.dtValue * Mn1V) *
                                  ((fbV - faV) * hn1V - (c_beta - c_alpha) * facV * hn1V +
                                   W * fbarriernV + nDependentMisfitAC1 + heterMechAC1));
  scalargradType  eqx_n1 = (constV(-userInputs.dtValue * Mn1V) * Knx1);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(0, eq_c);
  variable_list.set_scalar_gradient_term_RHS(0, eqx_c);

  variable_list.set_scalar_value_term_RHS(2, eq_n1);
  variable_list.set_scalar_gradient_term_RHS(2, eqx_n1);
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time
// independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are
// not explicit time-dependent equations. It takes "variable_list" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is
// given by "q_point_loc". The function outputs two terms to variable_list --
// one proportional to the test function and one proportional to the gradient of
// the test function. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const
{
  // --- Getting the values and derivatives of the model variables ---

  // The concentration and its derivatives (names here should match those in the
  // macros above)
  scalarvalueType c = variable_list.get_scalar_value(0);

  // The first order parameter and its derivatives (names here should match
  // those in the macros above)
  scalarvalueType n1 = variable_list.get_scalar_value(2);

  // The derivative of the displacement vector (names here should match those in
  // the macros above)
  vectorgradType ux = variable_list.get_vector_gradient(3);

  // --- Setting the expressions for the terms in the governing equations ---

  scalarvalueType h1V  = (3.0 * n1 * n1 - 2.0 * n1 * n1 * n1);
  scalarvalueType hn1V = (6.0 * n1 - 6.0 * n1 * n1);

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

  // This double-well function can be used to tune the interfacial energy
  // scalarvalueType fbarrierV = (n1*n1-2.0*n1*n1*n1+n1*n1*n1*n1);
  // scalarvalueType fbarriernV = (2.0*n1-6.0*n1*n1+4.0*n1*n1*n1);

  // Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
  scalarvalueType cacV, canV, cbnV, cbcV, cbcnV;

  cacV = fbccV / ((constV(1.0) - h1V) * fbccV + h1V * faccV);
  canV = hn1V * (c_alpha - c_beta) * cacV;

  cbcV  = faccV / ((constV(1.0) - h1V) * fbccV + h1V * faccV);
  cbnV  = hn1V * (c_alpha - c_beta) * cbcV;
  cbcnV = (faccV * (fbccV - faccV) * hn1V) /
          (((1.0 - h1V) * fbccV + h1V * faccV) *
           ((1.0 - h1V) * fbccV +
            h1V * faccV)); // Note: this is only true if faV and fbV are quadratic

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>> sfts1, sfts1c, sfts1cc, sfts1n,
    sfts1cn;

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
  dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E2[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]) - (sfts1[i][j] * h1V);
        }
    }

  // compute stress
  // S=C*(E-E0)
  //  Compute stress tensor (which is equal to the residual, Rux)
  dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

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

  vectorgradType eqx_u;

  // Fill residual corresponding to mechanics
  // R=-C*(E-E0)

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          eqx_u[i][j] = -S[i][j];
        }
    }

  scalarvalueType mu_c = constV(0.0);
  mu_c += facV * cacV * (constV(1.0) - h1V) + fbcV * cbcV * h1V;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          mu_c -= S[i][j] * (sfts1c[i][j] * h1V);
        }
    }

  scalarvalueType eq_mu = (mu_c);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(1, eq_mu);
  variable_list.set_vector_gradient_term_RHS(3, eqx_u);
}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. The
// (x,y,z) location of that quadrature point is given by "q_point_loc". The
// function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function -- for the
// left-hand-side of the equation. The index for each variable in this list
// corresponds to the index given at the top of this file. If there are multiple
// elliptic equations, conditional statements should be sed to ensure that the
// correct residual is being submitted. The index of the field being solved can
// be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const
{
  // --- Getting the values and derivatives of the model variables ---

  // n1
  scalarvalueType n1 = variable_list.get_scalar_value(2);

  // --- Setting the expressions for the terms in the governing equations ---

  vectorgradType eqx_Du;

  // Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the
  // dealii "symmetrize" function
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>> E;
  E = symmetrize(variable_list.get_change_in_vector_gradient(3));

  // Compute stress tensor (which is equal to the residual, Rux)
  if (n_dependent_stiffness == true)
    {
      scalarvalueType h1V = (3.0 * n1 * n1 - 2.0 * n1 * n1 * n1);
      dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double>> CIJ_combined;
      CIJ_combined = CIJ_Mg * (constV(1.0) - h1V);
      CIJ_combined += CIJ_Beta * (h1V);

      computeStress<dim>(CIJ_combined, E, eqx_Du);
    }
  else
    {
      computeStress<dim>(CIJ_Mg, E, eqx_Du);
    }

  // --- Submitting the terms for the governing equations ---

  variable_list.set_vector_gradient_term_LHS(3, eqx_Du);
}
