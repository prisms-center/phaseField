// List of residual equations for the coupled Allen-Cahn example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void
variableAttributeLoader::loadVariableAttributes()
{
  // Variable 1
  set_variable_name(0, "psi");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, IMPLICIT_TIME_DEPENDENT);

  set_dependencies_value_residual_term_RHS(0, "psi");
  set_dependencies_gradient_residual_term_RHS(0, "grad(psi)");
  set_dependencies_value_residual_term_LHS(0, "psi, change(psi)");
  set_dependencies_gradient_residual_term_LHS(0, "grad(change(psi))");
}

// =================================================================================
// residualRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of
// that quadrature point is given by "q_point_loc". The function outputs
// residuals to variable_list. The index for each variable in this list
// corresponds to the index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::residualExplicitRHS(
  variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const
{}

template <int dim, int degree>
void
customPDE<dim, degree>::residualNonexplicitRHS(
  variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const
{
  // The order parameter and its derivatives
  scalarvalueType psi  = variable_list.get_scalar_value(0);
  scalargradType  psix = variable_list.get_scalar_gradient(0);

  scalarvalueType psi_old = variable_list.get_old_scalar_value(0); // not implemented yet!

  // Parameters in the residual equations and expressions for the residual
  // equations can be set here.
  scalarvalueType fnV = (4.0 * psi * (psi - 1.0) * (psi - 0.5));

  // Residuals for the equation to evolve the order parameter
  scalarvalueType rpsi  = psi_old - psi - (constV(userInputs.dtValue * MnV) * fnV);
  scalargradType  rpsix = -(constV(userInputs.dtValue * KnV * MnV) * psix);

  // Submit the residuals
  variable_list.set_scalar_value_residual_term(0, rpsi);
  variable_list.set_scalar_gradient_residual_term(0, rpsix);
}

// =================================================================================
// residualLHS (needed only if at least one equation is elliptic)
// =================================================================================
// This function calculates the residual equations for the iterative solver for
// elliptic equations.for each variable. It takes "variable_list" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is
// given by "q_point_loc". The function outputs residual terms to
// "variable_list" for the left-hand-side of the residual equation for the
// iterative solver. The index for each variable in this list corresponds to the
// index given at the top of this file. If there are multiple elliptic
// equations, conditional statements should be used to ensure that the correct
// residual is being submitted. The index of the field being solved can be
// accessed by "this->currentFieldIndex".

template <int dim, int degree>
void
customPDE<dim, degree>::residualLHS(
  variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const
{
  // The order parameter and its derivatives

  scalarvalueType psi = variable_list.get_scalar_value(0);

  scalarvalueType Dpsi  = variable_list.get_change_in_scalar_value(0);
  scalargradType  Dpsix = variable_list.get_change_in_scalar_gradient(0);

  scalarvalueType rpsi =
    Dpsi + constV(userInputs.dtValue * MnV) *
             (12.0 * psi * psi * Dpsi - 12.0 * psi * Dpsi + 2.0 * Dpsi);
  scalargradType rpsix = constV(userInputs.dtValue * MnV * KnV) * Dpsix;

  variable_list.set_scalar_value_residual_term_LHS(0, rpsi);
  variable_list.set_scalar_gradient_residual_term_LHS(0, rpsix);
}
