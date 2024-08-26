// List of variables and residual equations for the Precipitate Evolution
// example application

// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 5

// The names of the variables, whether they are scalars or vectors and whether
// the governing eqn for the variable is parabolic or elliptic
#define variable_name          \
  {                            \
    "c", "n1", "n2", "n3", "u" \
  }
#define variable_type                                \
  {                                                  \
    "SCALAR", "SCALAR", "SCALAR", "SCALAR", "VECTOR" \
  }
#define variable_eq_type                                           \
  {                                                                \
    "PARABOLIC", "PARABOLIC", "PARABOLIC", "PARABOLIC", "ELLIPTIC" \
  }

// Flags for whether the value, gradient, and Hessian are needed in the residual
// eqns
#define need_val                  \
  {                               \
    true, true, true, true, false \
  }
#define need_grad                \
  {                              \
    true, true, true, true, true \
  }
#define need_hess                     \
  {                                   \
    false, false, false, false, false \
  }

// Flags for whether the residual equation has a term multiplied by the test
// function (need_val_residual) and/or the gradient of the test function
// (need_grad_residual)
#define need_val_residual         \
  {                               \
    true, true, true, true, false \
  }
#define need_grad_residual       \
  {                              \
    true, true, true, true, true \
  }

// Flags for whether the value, gradient, and Hessian are needed in the residual
// eqn for the left-hand-side of the iterative solver for elliptic equations
#define need_val_LHS               \
  {                                \
    false, true, true, true, false \
  }
#define need_grad_LHS                \
  {                                  \
    false, false, false, false, true \
  }
#define need_hess_LHS                 \
  {                                   \
    false, false, false, false, false \
  }

// Flags for whether the residual equation for the left-hand-side of the
// iterative solver for elliptic equations has a term multiplied by the test
// function (need_val_residual) and/or the gradient of the test function
// (need_grad_residual)
#define need_val_residual_LHS         \
  {                                   \
    false, false, false, false, false \
  }
#define need_grad_residual_LHS       \
  {                                  \
    false, false, false, false, true \
  }

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual
// equations can be set here. For simple cases, the entire residual equation can
// be written here. For more complex cases with loops or conditional statements,
// residual equations (or parts of residual equations) can be written below in
// "residualRHS".

// Cahn-Hilliard mobility
#define McV 1.0

// Allen-Cahn mobilities
#define Mn1V 50.0
#define Mn2V 50.0
#define Mn3V 50.0

// Gradient energy coefficients
double Kn1[3][3] = {
  {1.0, 0,   0  },
  {0,   0.5, 0  },
  {0,   0,   1.0}
};
double Kn2[3][3] = {
  {1.0, 0,   0  },
  {0,   0.5, 0  },
  {0,   0,   1.0}
};
double Kn3[3][3] = {
  {1.0, 0,   0  },
  {0,   0.5, 0  },
  {0,   0,   1.0}
};

// define energy barrier coefficient (used to tune the interfacial energy)
#define W -1.0

// Define mechanical properties
#define n_dependent_stiffness false
// Mechanical symmetry of the material and stiffness parameters
// If n_dependent_stiffness == false the first entry is used for all phases
#define MaterialModels \
  {                    \
    {                  \
      "ANISOTROPIC"    \
    }                  \
  }
#define MaterialConstants                                                              \
  {                                                                                    \
    {                                                                                  \
      31.3, 31.3, 32.45, 6.65, 6.65, 9.15, 13.0, 10.45, 0, 0, 0, 10.45, 0, 0, 0, 0, 0, \
        0, 0, 0, 0                                                                     \
    }                                                                                  \
  }

// Stress-free transformation strains
// Linear fits for the stress-free transformation strains in for sfts =
// sfts_linear * c + sfts_const
double sfts_linear1[3][3] = {
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};
double sfts_const1[3][3] = {
  {0.1305, 0,       0     },
  {0,      -0.0152, 0     },
  {0,      0,       -0.014}
}; // Mg-Nd beta-prime

double sfts_linear2[3][3] = {
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};
double sfts_const2[3][3] = {
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

double sfts_linear3[3][3] = {
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};
double sfts_const3[3][3] = {
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

// Free energy expressions
#define faV (24.7939 * c * c - 1.6752 * c - 1.9453e-06)
#define facV (49.5878 * c - 1.6752)
#define faccV (49.5878)
#define fbV (37.9316 * c * c - 10.7373 * c + 0.5401)
#define fbcV (75.8633 * c - 10.7373)
#define fbccV (75.8633)
#define h1V (3.0 * n1 * n1 - 2.0 * n1 * n1 * n1)
#define h2V (3.0 * n2 * n2 - 2.0 * n2 * n2 * n2)
#define h3V (3.0 * n3 * n3 - 2.0 * n3 * n3 * n3)
#define hn1V (6.0 * n1 - 6.0 * n1 * n1)
#define hn2V (6.0 * n2 - 6.0 * n2 * n2)
#define hn3V (6.0 * n3 - 6.0 * n3 * n3)

// This double-well function can be used to tune the interfacial energy
#define fbarrierV (n1 * n1 - 2.0 * n1 * n1 * n1 + n1 * n1 * n1 * n1)
#define fbarriernV (2.0 * n1 - 6.0 * n1 * n1 + 4.0 * n1 * n1 * n1)

// Residuals
#define rcV (c)
#define rcxTemp                                                         \
  (cx * ((1.0 - h1V - h2V - h3V) * faccV + (h1V + h2V + h3V) * fbccV) + \
   n1x * ((fbcV - facV) * hn1V) + n2x * ((fbcV - facV) * hn2V) +        \
   n3x * ((fbcV - facV) * hn3V))
#define rcxV (constV(-timeStep * McV) * rcxTemp)

#define rn1V                      \
  (n1 - constV(timeStep * Mn1V) * \
          ((fbV - faV) * hn1V + W * fbarriernV + nDependentMisfitAC1))
#define rn2V (n2 - constV(timeStep * Mn2V) * ((fbV - faV) * hn2V))
#define rn3V (n3 - constV(timeStep * Mn3V) * ((fbV - faV) * hn3V))
#define rn1xV (constV(-timeStep * Mn1V) * Knx1)
#define rn2xV (constV(-timeStep * Mn2V) * Knx2)
#define rn3xV (constV(-timeStep * Mn3V) * Knx3)

// =================================================================================
// residualRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "modelVariablesList" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. The
// (x,y,z) location of that quadrature point is given by "q_point_loc". The
// function outputs "modelResidualsList", a list of the value and gradient terms
// of the residual for each residual equation. The index for each variable in
// these lists corresponds to the order it is defined at the top of this file
// (starting at 0).
template <int dim>
void
generalizedProblem<dim>::residualRHS(
  const std::vector<modelVariable<dim>>              &modelVariablesList,
  std::vector<modelResidual<dim>>                    &modelResidualsList,
  dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const
{
  // The concentration and its derivatives (names here should match those in the
  // macros above)
  scalarvalueType c  = modelVariablesList[0].scalarValue;
  scalargradType  cx = modelVariablesList[0].scalarGrad;

  // The first order parameter and its derivatives (names here should match
  // those in the macros above)
  scalarvalueType n1  = modelVariablesList[1].scalarValue;
  scalargradType  n1x = modelVariablesList[1].scalarGrad;

  // The second order parameter and its derivatives (names here should match
  // those in the macros above)
  scalarvalueType n2  = modelVariablesList[2].scalarValue;
  scalargradType  n2x = modelVariablesList[2].scalarGrad;

  // The third order parameter and its derivatives (names here should match
  // those in the macros above)
  scalarvalueType n3  = modelVariablesList[3].scalarValue;
  scalargradType  n3x = modelVariablesList[3].scalarGrad;

  // The derivative of the displacement vector (names here should match those in
  // the macros above)
  vectorgradType ux = modelVariablesList[4].vectorGrad;
  vectorgradType ruxV;

  vectorhessType uxx;

  if (c_dependent_misfit == true)
    {
      uxx = modelVariablesList[4].vectorHess;
    }

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double>> sfts1, sfts1c, sfts1cc,
    sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

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
                CIJ_list[0][i][j] * (constV(1.0) - sum_hV) + CIJ_list[1][i][j] * sum_hV;
            }
        }
      computeStress<dim>(CIJ_combined, E2, S);
    }
  else
    {
      computeStress<dim>(CIJ_list[0], E2, S);
    }

  // Fill residual corresponding to mechanics
  // R=-C*(E-E0)

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          ruxV[i][j] = -S[i][j];
        }
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
      computeStress<dim>(CIJ_list[1] - CIJ_list[0], E2, S2);
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
          computeStress<dim>(CIJ_list[0], E3, S3);
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
                      grad_mu_el[k] += -S2[i][j] * (sfts1c[i][j] * hn1V * n1x[k] +
                                                    sfts2c[i][j] * hn2V * n2x[k] +
                                                    sfts3c[i][j] * hn3V * n3x[k]);
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

  modelResidualsList[0].scalarValueResidual = rcV;
  modelResidualsList[0].scalarGradResidual  = rcxV;

  modelResidualsList[1].scalarValueResidual = rn1V;
  modelResidualsList[1].scalarGradResidual  = rn1xV;

  modelResidualsList[2].scalarValueResidual = rn2V;
  modelResidualsList[2].scalarGradResidual  = rn2xV;

  modelResidualsList[3].scalarValueResidual = rn3V;
  modelResidualsList[3].scalarGradResidual  = rn3xV;

  modelResidualsList[4].vectorGradResidual = ruxV;
}

// =================================================================================
// residualLHS (needed only if at least one equation is elliptic)
// =================================================================================
// This function calculates the residual equations for the iterative solver for
// elliptic equations.for each variable. It takes "modelVariablesList" as an
// input, which is a list of the value and derivatives of each of the variables
// at a specific quadrature point. The (x,y,z) location of that quadrature point
// is given by "q_point_loc". The function outputs "modelRes", the value and
// gradient terms of for the left-hand-side of the residual equation for the
// iterative solver. The index for each variable in these lists corresponds to
// the order it is defined at the top of this file (starting at 0), not counting
// variables that have "need_val_LHS", "need_grad_LHS", and "need_hess_LHS" all
// set to "false". If there are multiple elliptic equations, conditional
// statements should be used to ensure that the correct residual is being
// submitted. The index of the field being solved can be accessed by
// "this->currentFieldIndex".
template <int dim>
void
generalizedProblem<dim>::residualLHS(
  const std::vector<modelVariable<dim>>              &modelVariablesList,
  modelResidual<dim>                                 &modelRes,
  dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const
{
  // n1
  scalarvalueType n1 = modelVariablesList[0].scalarValue;

  // n2
  scalarvalueType n2 = modelVariablesList[1].scalarValue;

  // n3
  scalarvalueType n3 = modelVariablesList[2].scalarValue;

  // u
  vectorgradType ux = modelVariablesList[3].vectorGrad;
  vectorgradType ruxV;

  // Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the
  // dealii "symmetrize" function
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>> E;
  E = symmetrize(ux);

  // Compute stress tensor (which is equal to the residual, Rux)
  if (n_dependent_stiffness == true)
    {
      dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double>> CIJ_combined;
      CIJ_combined =
        CIJ_list[0] * (constV(1.0) - h1V - h2V - h3V) + CIJ_list[1] * (h1V + h2V + h3V);

      computeStress<dim>(CIJ_combined, E, ruxV);
    }
  else
    {
      computeStress<dim>(CIJ_list[0], E, ruxV);
    }

  modelRes.vectorGradResidual = ruxV;
}

// =================================================================================
// energyDensity (needed only if calcEnergy == true)
// =================================================================================
// This function integrates the free energy density across the computational
// domain. It takes "modelVariablesList" as an input, which is a list of the
// value and derivatives of each of the variables at a specific quadrature
// point. It also takes the mapped quadrature weight, "JxW_value", as an input.
// The (x,y,z) location of the quadrature point is given by "q_point_loc". The
// weighted value of the energy density is added to "energy" variable and the
// components of the energy density are added to the "energy_components"
// variable (index 0: chemical energy, index 1: gradient energy, index 2:
// elastic energy).
template <int dim>
void
generalizedProblem<dim>::energyDensity(
  const std::vector<modelVariable<dim>>              &modelVarList,
  const dealii::VectorizedArray<double>              &JxW_value,
  dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc)
{
  scalarvalueType total_energy_density = constV(0.0);

  // c
  scalarvalueType c  = modelVarList[0].scalarValue;
  scalargradType  cx = modelVarList[0].scalarGrad;

  // n1
  scalarvalueType n1  = modelVarList[1].scalarValue;
  scalargradType  n1x = modelVarList[1].scalarGrad;

  // n2
  scalarvalueType n2  = modelVarList[2].scalarValue;
  scalargradType  n2x = modelVarList[2].scalarGrad;

  // n3
  scalarvalueType n3  = modelVarList[3].scalarValue;
  scalargradType  n3x = modelVarList[3].scalarGrad;

  // u
  vectorgradType ux = modelVarList[4].vectorGrad;

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
#if num_sop > 1
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += constV(0.5 * Kn2[i][j]) * n2x[i] * n2x[j];
        }
    }
#endif
#if num_sop > 2
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += constV(0.5 * Kn3[i][j]) * n3x[i] * n3x[j];
        }
    }
#endif

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double>> sfts1, sfts1c, sfts1cc,
    sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

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
          // E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V +
          // sfts2[i][j]*h2V + sfts3[i][j]*h3V);
          E2[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]) -
                     (sfts1[i][j] * h1V + sfts2[i][j] * h2V + sfts3[i][j] * h3V);
        }
    }

  // compute stress
  // S=C*(E-E0)
  dealii::VectorizedArray<double> CIJ_combined[2 * dim - 1 + dim / 3]
                                              [2 * dim - 1 + dim / 3];

  if (n_dependent_stiffness == true)
    {
      dealii::VectorizedArray<double> sum_hV;
      sum_hV = h1V + h2V + h3V;
      for (unsigned int i = 0; i < 2 * dim - 1 + dim / 3; i++)
        {
          for (unsigned int j = 0; j < 2 * dim - 1 + dim / 3; j++)
            {
              CIJ_combined[i][j] =
                CIJ_list[0][i][j] * (constV(1.0) - sum_hV) + CIJ_list[1][i][j] * sum_hV;
            }
        }
      computeStress<dim>(CIJ_combined, E2, S);
    }
  else
    {
      computeStress<dim>(CIJ_list[0], E2, S);
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

  assembler_lock.acquire();
  for (unsigned i = 0; i < c.size(); i++)
    {
      // For some reason, some of the values in this loop
      if (c[i] > 1.0e-10)
        {
          this->energy += total_energy_density[i] * JxW_value[i];
          this->energy_components[0] += f_chem[i] * JxW_value[i];
          this->energy_components[1] += f_grad[i] * JxW_value[i];
          this->energy_components[2] += f_el[i] * JxW_value[i];
        }
    }
  assembler_lock.release();
}
