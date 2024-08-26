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
  set_dependencies_gradient_term_RHS(
    0,
    "c, grad(c), n1, grad(n1), n2, grad(n2), n3, grad(n3), grad(u), hess(u)");

  // Variable 1
  set_variable_name(1, "n1");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(1, "c, n1, n2, n3, grad(u)");
  set_dependencies_gradient_term_RHS(1, "grad(n1)");

  // Variable 2
  set_variable_name(2, "n2");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(2, "c, n1, n2, n3, grad(u)");
  set_dependencies_gradient_term_RHS(2, "grad(n2)");

  // Variable 3
  set_variable_name(3, "n3");
  set_variable_type(3, SCALAR);
  set_variable_equation_type(3, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(3, "c, n1, n2, n3, grad(u)");
  set_dependencies_gradient_term_RHS(3, "grad(n3)");

  // Variable 2
  set_variable_name(4, "u");
  set_variable_type(4, VECTOR);
  set_variable_equation_type(4, TIME_INDEPENDENT);

  set_dependencies_value_term_RHS(4, "");
  set_dependencies_gradient_term_RHS(4, "c, n1, n2, n3, grad(u)");
  set_dependencies_value_term_LHS(4, "");
  set_dependencies_gradient_term_LHS(4, "n1, n2, n3, grad(change(u))");
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

  // The concentration and its derivatives
  scalarvalueType c  = variable_list.get_scalar_value(0);
  scalargradType  cx = variable_list.get_scalar_gradient(0);

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

  // --- Setting the expressions for the terms in the governing equations ---

  vectorhessType uxx;

  bool c_dependent_misfit = false;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          if (std::abs(sfts_linear1[i][j]) > 1.0e-12)
            {
              c_dependent_misfit = true;
            }
        }
    }

  if (c_dependent_misfit == true)
    {
      uxx = variable_list.get_vector_hessian(4);
    }

  // Free energy expressions and interpolation functions
  scalarvalueType faV  = (A0 + A1 * c + A2 * c * c + A3 * c * c * c + A4 * c * c * c * c);
  scalarvalueType facV = (A1 + 2.0 * A2 * c + 3.0 * A3 * c * c + 4.0 * A4 * c * c * c);
  scalarvalueType faccV = (2.0 * A2 + 6.0 * A3 * c + 12.0 * A4 * c * c);
  scalarvalueType fbV   = (B2 * c * c + B1 * c + B0);
  scalarvalueType fbcV  = (2.0 * B2 * c + B1);
  scalarvalueType fbccV = constV(2.0 * B2);
  scalarvalueType h1V =
    (10.0 * n1 * n1 * n1 - 15.0 * n1 * n1 * n1 * n1 + 6.0 * n1 * n1 * n1 * n1 * n1);
  scalarvalueType h2V =
    (10.0 * n2 * n2 * n2 - 15.0 * n2 * n2 * n2 * n2 + 6.0 * n2 * n2 * n2 * n2 * n2);
  scalarvalueType h3V =
    (10.0 * n3 * n3 * n3 - 15.0 * n3 * n3 * n3 * n3 + 6.0 * n3 * n3 * n3 * n3 * n3);
  scalarvalueType hn1V =
    (30.0 * n1 * n1 - 60.0 * n1 * n1 * n1 + 30.0 * n1 * n1 * n1 * n1);
  scalarvalueType hn2V =
    (30.0 * n2 * n2 - 60.0 * n2 * n2 * n2 + 30.0 * n2 * n2 * n2 * n2);
  scalarvalueType hn3V =
    (30.0 * n3 * n3 - 60.0 * n3 * n3 * n3 + 30.0 * n3 * n3 * n3 * n3);

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>> sfts1, sfts1c, sfts1cc, sfts2,
    sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

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
  dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

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
  //  Compute stress tensor (which is equal to the residual, Rux)
  dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

  if (n_dependent_stiffness == true)
    {
      dealii::VectorizedArray<double> sum_hV;
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

  // Compute one of the stress terms in the order parameter chemical potential,
  // nDependentMisfitACp = C*(E-E0)*(E0_p*Hn)
  dealii::VectorizedArray<double> nDependentMisfitAC1 = constV(0.0);
  dealii::VectorizedArray<double> nDependentMisfitAC2 = constV(0.0);
  dealii::VectorizedArray<double> nDependentMisfitAC3 = constV(0.0);

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          nDependentMisfitAC1 += S[i][j] * (sfts1[i][j]);
          nDependentMisfitAC2 += S[i][j] * (sfts2[i][j]);
          nDependentMisfitAC3 += S[i][j] * (sfts3[i][j]);
        }
    }

  nDependentMisfitAC1 *= -hn1V;
  nDependentMisfitAC2 *= -hn2V;
  nDependentMisfitAC3 *= -hn3V;

  // Compute the other stress term in the order parameter chemical potential,
  // heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
  dealii::VectorizedArray<double> heterMechAC1 = constV(0.0);
  dealii::VectorizedArray<double> heterMechAC2 = constV(0.0);
  dealii::VectorizedArray<double> heterMechAC3 = constV(0.0);
  dealii::VectorizedArray<double> S2[dim][dim];

  if (n_dependent_stiffness == true)
    {
      // computeStress<dim>(CIJ_diff, E2, S2);
      computeStress<dim>(CIJ_Beta - CIJ_Mg, E2, S2);
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              heterMechAC1 += S2[i][j] * E2[i][j];
            }
        }
      // Aside from HnpV, heterMechAC1, heterMechAC2, and heterMechAC3 are equal
      heterMechAC2 = 0.5 * hn2V * heterMechAC1;
      heterMechAC3 = 0.5 * hn3V * heterMechAC1;

      heterMechAC1 = 0.5 * hn1V * heterMechAC1;
    }

  // compute the stress term in the gradient of the concentration chemical
  // potential, grad_mu_el = [C*(E-E0)*E0c]x, must be a vector with length dim
  scalargradType grad_mu_el;

  if (c_dependent_misfit == true)
    {
      dealii::VectorizedArray<double> E3[dim][dim], S3[dim][dim];

      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              E3[i][j] = -(sfts1c[i][j] * h1V + sfts2c[i][j] * h2V + sfts3c[i][j] * h3V);
            }
        }

      if (n_dependent_stiffness == true)
        {
          computeStress<dim>(CIJ_combined, E3, S3);
        }
      else
        {
          computeStress<dim>(CIJ_Mg, E3, S3);
        }

      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              for (unsigned int k = 0; k < dim; k++)
                {
                  grad_mu_el[k] +=
                    S3[i][j] *
                    (constV(0.5) * (uxx[i][j][k] + uxx[j][i][k]) + E3[i][j] * cx[k] -
                     (sfts1[i][j] * hn1V * n1x[k] + sfts2[i][j] * hn2V * n2x[k] +
                      sfts3[i][j] * hn3V * n3x[k]));

                  grad_mu_el[k] +=
                    -S[i][j] *
                    (sfts1c[i][j] * hn1V * n1x[k] + sfts2c[i][j] * hn2V * n2x[k] +
                     sfts3c[i][j] * hn3V * n3x[k] +
                     (sfts1cc[i][j] * h1V + sfts2cc[i][j] * h2V + sfts3cc[i][j] * h3V) *
                       cx[k]);

                  if (n_dependent_stiffness == true)
                    {
                      grad_mu_el[k] += S2[i][j] * E3[i][j] *
                                       (hn1V * n1x[k] + hn2V * n2x[k] + hn3V * n3x[k]);
                    }
                }
            }
        }
    }

  // compute K*nx
  scalargradType Knx1, Knx2, Knx3;
  for (unsigned int a = 0; a < dim; a++)
    {
      Knx1[a] = 0.0;
      Knx2[a] = 0.0;
      Knx3[a] = 0.0;
      for (unsigned int b = 0; b < dim; b++)
        {
          Knx1[a] += constV(Kn1[a][b]) * n1x[b];
          Knx2[a] += constV(Kn2[a][b]) * n2x[b];
          Knx3[a] += constV(Kn3[a][b]) * n3x[b];
        }
    }

  // The terms in the govering equations
  scalarvalueType eq_c = (c);
  scalargradType  eqx_c_temp =
    (cx * ((1.0 - h1V - h2V - h3V) * faccV + (h1V + h2V + h3V) * fbccV) +
     n1x * ((fbcV - facV) * hn1V) + n2x * ((fbcV - facV) * hn2V) +
     n3x * ((fbcV - facV) * hn3V) + grad_mu_el);
  scalargradType eqx_c = (constV(-userInputs.dtValue * McV) * eqx_c_temp);

  scalarvalueType eq_n1 =
    (n1 - constV(userInputs.dtValue * Mn1V) *
            ((fbV - faV) * hn1V + nDependentMisfitAC1 + heterMechAC1));
  scalarvalueType eq_n2 =
    (n2 - constV(userInputs.dtValue * Mn2V) *
            ((fbV - faV) * hn2V + nDependentMisfitAC2 + heterMechAC2));
  scalarvalueType eq_n3 =
    (n3 - constV(userInputs.dtValue * Mn3V) *
            ((fbV - faV) * hn3V + nDependentMisfitAC3 + heterMechAC3));
  scalargradType eqx_n1 = (constV(-userInputs.dtValue * Mn1V) * Knx1);
  scalargradType eqx_n2 = (constV(-userInputs.dtValue * Mn2V) * Knx2);
  scalargradType eqx_n3 = (constV(-userInputs.dtValue * Mn3V) * Knx3);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(0, eq_c);
  variable_list.set_scalar_gradient_term_RHS(0, eqx_c);

  variable_list.set_scalar_value_term_RHS(1, eq_n1);
  variable_list.set_scalar_gradient_term_RHS(1, eqx_n1);

  variable_list.set_scalar_value_term_RHS(2, eq_n2);
  variable_list.set_scalar_gradient_term_RHS(2, eqx_n2);

  variable_list.set_scalar_value_term_RHS(3, eq_n3);
  variable_list.set_scalar_gradient_term_RHS(3, eqx_n3);
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

  // The concentration and its derivatives
  scalarvalueType c = variable_list.get_scalar_value(0);

  // The first order parameter and its derivatives
  scalarvalueType n1 = variable_list.get_scalar_value(1);

  // The second order parameter and its derivatives
  scalarvalueType n2 = variable_list.get_scalar_value(2);

  // The third order parameter and its derivatives
  scalarvalueType n3 = variable_list.get_scalar_value(3);

  // The derivative of the displacement vector
  vectorgradType ux = variable_list.get_vector_gradient(4);

  // --- Setting the expressions for the terms in the governing equations ---

  // Interpolation functions
  scalarvalueType h1V =
    (10.0 * n1 * n1 * n1 - 15.0 * n1 * n1 * n1 * n1 + 6.0 * n1 * n1 * n1 * n1 * n1);
  scalarvalueType h2V =
    (10.0 * n2 * n2 * n2 - 15.0 * n2 * n2 * n2 * n2 + 6.0 * n2 * n2 * n2 * n2 * n2);
  scalarvalueType h3V =
    (10.0 * n3 * n3 * n3 - 15.0 * n3 * n3 * n3 * n3 + 6.0 * n3 * n3 * n3 * n3 * n3);

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>> sfts1, sfts2, sfts3;

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts1[i][j] = constV(sfts_linear1[i][j]) * c + constV(sfts_const1[i][j]);
          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts2[i][j] = constV(sfts_linear2[i][j]) * c + constV(sfts_const2[i][j]);
          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts3[i][j] = constV(sfts_linear3[i][j]) * c + constV(sfts_const3[i][j]);
        }
    }

  // compute E2=(E-E0)
  dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

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
  //  Compute stress tensor (which is equal to the residual, Rux)
  dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

  if (n_dependent_stiffness == true)
    {
      dealii::VectorizedArray<double> sum_hV;
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

  vectorgradType eqx_u;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          eqx_u[i][j] = -S[i][j];
        }
    }

  // --- Submitting the terms for the governing equations ---

  variable_list.set_vector_gradient_term_RHS(4, eqx_u);
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
  scalarvalueType n1 = variable_list.get_scalar_value(1);

  // n2
  scalarvalueType n2 = variable_list.get_scalar_value(2);

  // n3
  scalarvalueType n3 = variable_list.get_scalar_value(3);

  // u
  vectorgradType Dux = variable_list.get_change_in_vector_gradient(4);

  // --- Setting the expressions for the terms in the governing equations ---

  vectorgradType eqx_Du;

  // Interpolation functions

  scalarvalueType h1V =
    (10.0 * n1 * n1 * n1 - 15.0 * n1 * n1 * n1 * n1 + 6.0 * n1 * n1 * n1 * n1 * n1);
  scalarvalueType h2V =
    (10.0 * n2 * n2 * n2 - 15.0 * n2 * n2 * n2 * n2 + 6.0 * n2 * n2 * n2 * n2 * n2);
  scalarvalueType h3V =
    (10.0 * n3 * n3 * n3 - 15.0 * n3 * n3 * n3 * n3 + 6.0 * n3 * n3 * n3 * n3 * n3);

  // Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the
  // dealii "symmetrize" function
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>> E;
  E = symmetrize(Dux);

  // Compute stress tensor (which is equal to the residual, Rux)
  if (n_dependent_stiffness == true)
    {
      dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double>> CIJ_combined;
      CIJ_combined =
        CIJ_Mg * (constV(1.0) - h1V - h2V - h3V) + CIJ_Beta * (h1V + h2V + h3V);

      computeStress<dim>(CIJ_combined, E, eqx_Du);
    }
  else
    {
      computeStress<dim>(CIJ_Mg, E, eqx_Du);
    }

  // --- Submitting the terms for the governing equations ---

  variable_list.set_vector_gradient_term_LHS(4, eqx_Du);
}
