// List of variables and residual equations for the mechanics example
// application

// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 4

// The names of the variables, whether they are scalars or vectors and whether
// the governing eqn for the variable is parabolic or elliptic
#define variable_name    \
  {                      \
    "u", "u2", "u3", "c" \
  }
#define variable_type                      \
  {                                        \
    "VECTOR", "VECTOR", "VECTOR", "SCALAR" \
  }
#define variable_eq_type                            \
  {                                                 \
    "ELLIPTIC", "PARABOLIC", "ELLIPTIC", "ELLIPTIC" \
  }

// Flags for whether the value, gradient, and Hessian are needed in the residual
// eqns
#define need_val               \
  {                            \
    false, false, false, false \
  }
#define need_grad          \
  {                        \
    true, true, true, true \
  }
#define need_hess              \
  {                            \
    false, false, false, false \
  }

// Flags for whether the residual equation has a term multiplied by the test
// function (need_val_residual) and/or the gradient of the test function
// (need_grad_residual)
#define need_val_residual      \
  {                            \
    false, false, false, false \
  }
#define need_grad_residual \
  {                        \
    true, true, true, true \
  }

// Flags for whether the value, gradient, and Hessian are needed in the residual
// eqn for the left-hand-side of the iterative solver for elliptic equations
#define need_val_LHS           \
  {                            \
    false, false, false, false \
  }
#define need_grad_LHS      \
  {                        \
    true, true, true, true \
  }
#define need_hess_LHS          \
  {                            \
    false, false, false, false \
  }

// Flags for whether the residual equation for the left-hand-side of the
// iterative solver for elliptic equations has a term multiplied by the test
// function (need_val_residual) and/or the gradient of the test function
// (need_grad_residual)
#define need_val_residual_LHS  \
  {                            \
    false, false, false, false \
  }
#define need_grad_residual_LHS \
  {                            \
    true, false, true, true    \
  }

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual
// equations can be set here. For simple cases, the entire residual equation can
// be written here. For more complex cases with loops or conditional statements,
// residual equations (or parts of residual equations) can be written below in
// "residualRHS".

// Define Mechanical properties
// Mechanical symmetry of the material and stiffness parameters
#define MaterialModels \
  {                    \
    "ISOTROPIC"        \
  }
#define MaterialConstants \
  {                       \
    {                     \
      2.0, 0.3            \
    }                     \
  }

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
  // u
  vectorgradType ux  = modelVariablesList[0].vectorGrad;
  vectorgradType u2x = modelVariablesList[1].vectorGrad;
  vectorgradType u3x = modelVariablesList[2].vectorGrad;
  scalargradType cx  = modelVariablesList[3].scalarGrad;
  vectorgradType Rux, Rux2;

  // compute strain tensor
  dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]);
        }
    }

  // compute strain tensor
  dealii::VectorizedArray<double> E2[dim][dim], S2[dim][dim];
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E2[i][j] = constV(0.5) * (u3x[i][j] + u3x[j][i]);
        }
    }

  // compute stress tensor
  computeStress<dim>(CIJ_list[0], E, S);

  computeStress<dim>(CIJ_list[0], E2, S2);

  // compute residual
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          Rux[i][j]  = -S[i][j];
          Rux2[i][j] = -S2[i][j];
        }
    }

  modelResidualsList[0].vectorGradResidual = Rux;
  modelResidualsList[1].vectorGradResidual = constV(0.0) * Rux;
  modelResidualsList[2].vectorGradResidual = Rux2;
  modelResidualsList[3].scalarGradResidual = -cx;
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
  const std::vector<modelVariable<dim>>              &modelVarList,
  modelResidual<dim>                                 &modelRes,
  dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const
{
  // u
  vectorgradType ux  = modelVarList[0].vectorGrad;
  vectorgradType u2x = modelVarList[1].vectorGrad;
  vectorgradType u3x = modelVarList[2].vectorGrad;
  scalargradType cx  = modelVarList[3].scalarGrad;
  vectorgradType Rux;

  dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
  if (this->currentFieldIndex == 0)
    {
      // compute strain tensor
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              E[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]);
            }
        }
    }
  else
    {
      // compute strain tensor
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              E[i][j] = constV(0.5) * (u3x[i][j] + u3x[j][i]);
            }
        }
    }

  // compute stress tensor
  computeStress<dim>(CIJ_list[0], E, S);

  // compute residual
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          Rux[i][j] = S[i][j];
        }
    }

  if (this->currentFieldIndex < 3)
    {
      modelRes.vectorGradResidual = Rux;
    }
  else
    {
      modelRes.scalarGradResidual = cx;
    }
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
{}
