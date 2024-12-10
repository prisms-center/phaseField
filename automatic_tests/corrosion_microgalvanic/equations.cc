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
#include <cmath>
#include <iostream>
using namespace std;

void
variableAttributeLoader::loadVariableAttributes()
{
  // Declaring all the model variables as follows:
  // domain parameter for the electrolyte - psi
  // order parameters for the anodic phase, nAnodic, and the cathodic phase,
  // nCathodic
  // chemical potentials for the anodic phase, muAnodic,
  // and the cathodic phase,
  // muCathodic
  // electrolyte potential - Phi
  // reaction current density - irxn
  // interpolation factor for the anodic phase, xiAnodic

  // Variable 0
  set_variable_name(0, "nAnodic");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "nAnodic, irxn, grad(psi), Phi, xiAnodic");
  set_dependencies_gradient_term_RHS(0, "nAnodic, grad(muAnodic), irxn, Phi, xiAnodic");

  // Variable 1
  set_variable_name(1, "muAnodic");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, AUXILIARY);

  set_dependencies_value_term_RHS(1, "nAnodic, nCathodic, psi");
  set_dependencies_gradient_term_RHS(1, "grad(nAnodic)");

  // Variable 2
  set_variable_name(2, "nCathodic");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(2, "nCathodic");
  set_dependencies_gradient_term_RHS(2, "nCathodic, grad(muCathodic), irxn, Phi");

  // Variable 3
  set_variable_name(3, "muCathodic");
  set_variable_type(3, SCALAR);
  set_variable_equation_type(3, AUXILIARY);

  set_dependencies_value_term_RHS(3, "nCathodic, nAnodic, psi");
  set_dependencies_gradient_term_RHS(3, "grad(nCathodic)");

  // Variable 4
  set_variable_name(4, "psi");
  set_variable_type(4, SCALAR);
  set_variable_equation_type(4, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(4, "psi, irxn, grad(psi), Phi, xiAnodic");
  set_dependencies_gradient_term_RHS(4, "psi, grad(mupsi), irxn, Phi, xiAnodic");

  // Variable 5
  set_variable_name(5, "mupsi");
  set_variable_type(5, SCALAR);
  set_variable_equation_type(5, AUXILIARY);

  set_dependencies_value_term_RHS(5, "nAnodic, nCathodic, psi");
  set_dependencies_gradient_term_RHS(5, "grad(psi)");

  // Variable 6
  set_variable_name(6, "Phi");
  set_variable_type(6, SCALAR);
  set_variable_equation_type(6, TIME_INDEPENDENT);

  set_dependencies_value_term_LHS(
    6,
    "nAnodic, nCathodic, grad(psi), change(Phi), Phi, xiAnodic");
  set_dependencies_gradient_term_LHS(6, "psi, grad(change(Phi))");
  set_dependencies_value_term_RHS(6, "grad(psi), irxn, xiAnodic");
  set_dependencies_gradient_term_RHS(6, "psi, grad(Phi)");

  // Variable 7
  set_variable_name(7, "irxn");
  set_variable_type(7, SCALAR);
  set_variable_equation_type(7, AUXILIARY);

  set_dependencies_value_term_RHS(7, "nCathodic, nAnodic, Phi, xiAnodic");
  set_dependencies_gradient_term_RHS(7, "");

  // Variable 8
  set_variable_name(8, "xiAnodic");
  set_variable_type(8, SCALAR);
  set_variable_equation_type(8, AUXILIARY);

  set_dependencies_value_term_RHS(8, "nAnodic, nCathodic");
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
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // --- Parameters in the explicit equations can be set here  ---

  // Timestep
  scalarvalueType delt = constV(userInputs.dtValue);

  // The order parameter of the anodic phase
  scalarvalueType nAnodic = variable_list.get_scalar_value(0);

  // The gradient of the chemical potential of the anodic phase
  scalargradType muAnodicx = variable_list.get_scalar_gradient(1);

  // The order parameter of the cathodic phase
  scalarvalueType nCathodic = variable_list.get_scalar_value(2);

  // The gradient of the chemical potential of the cathodic phase
  scalargradType muCathodicx = variable_list.get_scalar_gradient(3);

  // The domain parameter of the electrolyte phase and its derivatives
  scalarvalueType psi  = variable_list.get_scalar_value(4);
  scalargradType  psix = variable_list.get_scalar_gradient(4);

  // The derivative of the chemical potential of the electrolyte phase
  scalargradType mupsix = variable_list.get_scalar_gradient(5);

  // The electrolyte potential
  scalarvalueType Phi = variable_list.get_scalar_value(6);

  // The interpolation factor for the anodic phase
  scalarvalueType xiAnodic = variable_list.get_scalar_value(8);

  // Ensuring that the domain and order parameters remain within [0 1]

  // Electrolyte
  psi = min(psi, constV(1.0));
  psi = max(psi, constV(0.0));

  // Anodic phase
  nAnodic = min(nAnodic, constV(1.0));
  nAnodic = max(nAnodic, constV(0.0));

  // Cathodic phase
  nCathodic = min(nCathodic, constV(1.0));
  nCathodic = max(nCathodic, constV(0.0));

  // --- Calculation of terms needed in multiple expressions  ---

  // Magnifude of the derivative of psi
  scalarvalueType magpsix = constV(0.0);
  for (int i = 0; i < dim; i++)
    {
      magpsix += psix[i] * psix[i];
    }
  magpsix = sqrt(magpsix);

  // --- Calculation of irxn  ---

  // Overpotentials for the anodic and cathodic phases
  scalarvalueType etaAnodic = VsV - EcorrAnodic - Phi;

  // Calculation of anodic/cathodic current density
  scalarvalueType itafel  = i0Anodic * exp(etaAnodic / AAnodic);
  scalarvalueType iAnodic = itafel / (constV(1.0) + (itafel / iMax));

  // Time for initial equilibration of all the interfaces present in the
  // domain

  // Present simulation time
  double simTime = this->currentTime;

  // Turning on the corrosion-caused motion of the interface after
  // tStepsForV * timeStep time. The number 0.25 decides the width of the
  // smooth transition between vMin and vMax, as described below
  double arg = (simTime - tStepStartForV * userInputs.dtValue) / 0.25;

  // switchVal is zero before tStepsForV * timeStep and 1 after that time
  double switchVal = 0.5 * (1.0 + tanh(arg));

  // Normal component of the interface velocity, which switches from vMin
  // to vMax at time = tStepStartForV * timeStep

  scalarvalueType vMin = constV(0.0);
  scalarvalueType vMax = -constV(VMV / (zMV * FarC)) * xiAnodic * iAnodic;
  scalarvalueType v    = constV(1.0 - switchVal) * vMin + constV(switchVal) * vMax;

  // The mobility for all the phases
  scalarvalueType MnV = MconstV * psi * abs(xiAnodic * iAnodic);

  // --- Calculation of residual terms for nAnodic  ---
  scalarvalueType rnAnodic  = nAnodic + v * delt * magpsix;
  scalargradType  rnAnodicx = -MnV * delt * muAnodicx;

  // --- Calculation of residual terms for nCathodic  ---
  scalarvalueType rnCathodic  = nCathodic;
  scalargradType  rnCathodicx = -MnV * delt * muCathodicx;

  // --- Calculation of residual terms for psi ---
  scalarvalueType rpsi  = psi - v * delt * magpsix;
  scalargradType  rpsix = -MnV * delt * mupsix;

  // --- Submitting the terms for the governing equations ---
  // Residuals terms for the equation to evolve the order parameter
  variable_list.set_scalar_value_term_RHS(0, rnAnodic);
  variable_list.set_scalar_gradient_term_RHS(0, rnAnodicx);

  // Residuals terms for the equation to evolve the order parameter
  variable_list.set_scalar_value_term_RHS(2, rnCathodic);
  variable_list.set_scalar_gradient_term_RHS(2, rnCathodicx);

  // Residuals terms for the equation to evolve the domain parameter
  variable_list.set_scalar_value_term_RHS(4, rpsi);
  variable_list.set_scalar_gradient_term_RHS(4, rpsix);
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
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // --- Getting the values and derivatives of the model variables ---

  // The order parameter of the anodic phase and its gradient
  scalarvalueType nAnodic  = variable_list.get_scalar_value(0);
  scalargradType  nAnodicx = variable_list.get_scalar_gradient(0);

  // The order parameter of the cathodic phase and its gradient
  scalarvalueType nCathodic  = variable_list.get_scalar_value(2);
  scalargradType  nCathodicx = variable_list.get_scalar_gradient(2);

  // The domain parameter of the electrolyte phase and its derivatives
  scalarvalueType psi  = variable_list.get_scalar_value(4);
  scalargradType  psix = variable_list.get_scalar_gradient(4);

  // The electrolyte potential and its derivatives
  scalarvalueType Phi  = variable_list.get_scalar_value(6);
  scalargradType  Phix = variable_list.get_scalar_gradient(6);

  // The reaction current density
  scalarvalueType irxn     = variable_list.get_scalar_value(7);
  scalarvalueType xiAnodic = variable_list.get_scalar_value(8);

  // Ensuring domain and order parameters are in the range [0 1]
  psi = min(psi, constV(1.0));
  psi = max(psi, constV(0.0));

  nAnodic = min(nAnodic, constV(1.0));
  nAnodic = max(nAnodic, constV(0.0));

  nCathodic = min(nCathodic, constV(1.0));
  nCathodic = max(nCathodic, constV(0.0));

  // --- Calculation of terms needed in multiple expressions  ---
  // Magnitude of the gradient of the domain parameter

  scalarvalueType magpsix = constV(0.0);
  for (int i = 0; i < dim; i++)
    {
      magpsix = magpsix + psix[i] * psix[i];
    }
  magpsix = sqrt(magpsix);

  // --- Calculation of residual terms for muAnodic  ---

  // Derivative of bulk free energy with respect to nAnodic
  scalarvalueType rmuAnodic =
    (nAnodic * nAnodic * nAnodic - nAnodic + constV(2.0) * gamma * nAnodic * psi * psi +
     constV(2.0) * gamma * nAnodic * nCathodic * nCathodic);

  scalargradType rmuAnodicx = epssqV * nAnodicx;

  // --- Calculation of residual terms for muCathodic  ---

  // Derivative of bulk free energy with respect to nCathodic
  scalarvalueType rmuCathodic = (nCathodic * nCathodic * nCathodic - nCathodic +
                                 constV(2.0) * gamma * nCathodic * psi * psi +
                                 constV(2.0) * gamma * nCathodic * nAnodic * nAnodic);

  scalargradType rmuCathodicx = epssqV * nCathodicx;

  // --- Calculation of residual terms for mupsi  ---

  // Derivative of bulk free energy with respect to psi
  scalarvalueType rmupsi =
    (psi * psi * psi - psi + constV(2.0) * gamma * psi * nAnodic * nAnodic +
     constV(2.0) * gamma * psi * nCathodic * nCathodic);

  scalargradType rmupsix = epssqV * psix;

  // --- Calculation of irxn  ---

  // Overpotential

  scalarvalueType etaAnodic   = VsV - EcorrAnodic - Phi;
  scalarvalueType etaCathodic = VsV - EcorrCathodic - Phi;

  // Calculation of anodic/cathodic current density
  scalarvalueType itafel    = i0Anodic * exp(etaAnodic / AAnodic);
  scalarvalueType iAnodic   = itafel / (constV(1.0) + (itafel / iMax));
  scalarvalueType iCathodic = -i0Cathodic * exp(etaCathodic / ACathodic);

  // Calculation of the interpolation factor for the anodic phase
  // The factor for the cathodic phase is calculated  first, and then the
  // one for the anodic phase is calculated via xiAnodic = 1 - xiCathodic
  scalarvalueType xiCathodic = nCathodic / max(nAnodic + nCathodic, constV(lthresh));

  xiAnodic = constV(1.0) - xiCathodic;

  // combination of irxn
  scalarvalueType rirxn = xiAnodic * iAnodic + xiCathodic * iCathodic;

  // --- Calculation of residual terms for Phi  ---
  scalarvalueType rPhi  = -magpsix * irxn;
  scalargradType  rPhix = psi * kappa * Phix;

  // --- Submitting the terms for the governing equations ---
  // Residuals for the equation to calculate muAnodic
  variable_list.set_scalar_value_term_RHS(1, rmuAnodic);
  variable_list.set_scalar_gradient_term_RHS(1, rmuAnodicx);

  // Residuals for the equation to calculate muCathodic
  variable_list.set_scalar_value_term_RHS(3, rmuCathodic);
  variable_list.set_scalar_gradient_term_RHS(3, rmuCathodicx);

  // Residuals for the equation to calculate mupsi
  variable_list.set_scalar_value_term_RHS(5, rmupsi);
  variable_list.set_scalar_gradient_term_RHS(5, rmupsix);

  // Residuals for the equation to evolve the Potential
  variable_list.set_scalar_value_term_RHS(6, rPhi);
  variable_list.set_scalar_gradient_term_RHS(6, rPhix);

  // irxn for the equation to evolve irxn
  variable_list.set_scalar_value_term_RHS(7, rirxn);

  // Equation to solve for xiAnodic
  variable_list.set_scalar_value_term_RHS(8, xiAnodic);
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
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // The order parameter of the anodic phase
  scalarvalueType nAnodic = variable_list.get_scalar_value(0);

  // The order parameter of the cathodic phase
  scalarvalueType nCathodic = variable_list.get_scalar_value(2);

  // The domain parameter of the electrolyt phase and its derivatives
  scalarvalueType psi  = variable_list.get_scalar_value(4);
  scalargradType  psix = variable_list.get_scalar_gradient(4);

  // The electrolyte potential and its derivatives
  scalarvalueType Phi = variable_list.get_scalar_value(6);

  // The change in potential in the electrode and its derivatives
  scalarvalueType DPhi  = variable_list.get_change_in_scalar_value(6);
  scalargradType  DPhix = variable_list.get_change_in_scalar_gradient(6);

  // The interpolation factor for the anodic phase
  scalarvalueType xiAnodic = variable_list.get_scalar_value(8);

  // Ensuring that the domain and order parameters remain within the range
  // [0 1]
  psi = min(psi, constV(1.0));
  psi = max(psi, constV(0.0));

  nAnodic = min(nAnodic, constV(1.0));
  nAnodic = max(nAnodic, constV(0.0));

  nCathodic = min(nCathodic, constV(1.0));
  nCathodic = max(nCathodic, constV(0.0));

  // --- Calculation of terms needed in multiple expressions  ---

  // Magnitude of the gradient of the domain parameter
  scalarvalueType magpsix = constV(0.0);
  for (int i = 0; i < dim; i++)
    {
      magpsix = magpsix + psix[i] * psix[i];
    }
  magpsix = sqrt(magpsix);

  // --- Calculation of residual terms for DPhi (delta phi on the LHS) ---
  // --- Calculation of residual terms for irxnp (the derivative of
  // irxn wrt to Phi)   ---

  // Overpotential
  scalarvalueType etaAnodic   = VsV - EcorrAnodic - Phi;
  scalarvalueType etaCathodic = VsV - EcorrCathodic - Phi;

  // Calculation of anodic/cathodic current density
  scalarvalueType itafel    = i0Anodic * exp(etaAnodic / AAnodic);
  scalarvalueType iCathodic = -i0Cathodic * exp(etaCathodic / ACathodic);

  // The interpolation factor for the cathodic phase
  scalarvalueType xiCathodic = constV(1.0) - xiAnodic;

  // The derivative of irxn wrt to Phi
  scalarvalueType irxnp =
    -xiAnodic * (iMax / (iMax + itafel)) * (iMax / (iMax + itafel)) * (itafel / AAnodic) -
    xiCathodic * (iCathodic / ACathodic);

  scalarvalueType rDPhi  = magpsix * irxnp * DPhi;
  scalargradType  rDPhix = -psi * kappa * DPhix;

  // Residuals for the equation to evolve the potential(Phi)
  variable_list.set_scalar_value_term_LHS(6, rDPhi);
  variable_list.set_scalar_gradient_term_LHS(6, rDPhix);
}
