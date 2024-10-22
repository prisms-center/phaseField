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
  set_variable_name(0, "n");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "n, grad(psi), irxn");
  set_dependencies_gradient_term_RHS(0, "psi, grad(mu), irxn");

  // Variable 1
  set_variable_name(1, "mu");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, AUXILIARY);

  set_dependencies_value_term_RHS(1, "n, psi");
  set_dependencies_gradient_term_RHS(1, "grad(n)");

  // Variable 2
  set_variable_name(2, "psi");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(2, "psi, grad(psi), irxn");
  set_dependencies_gradient_term_RHS(2, "psi, grad(mupsi), irxn");

  // Variable 3
  set_variable_name(3, "mupsi");
  set_variable_type(3, SCALAR);
  set_variable_equation_type(3, AUXILIARY);

  set_dependencies_value_term_RHS(3, "n, psi");
  set_dependencies_gradient_term_RHS(3, "grad(psi)");

  // Variable 4
  set_variable_name(4, "cM");
  set_variable_type(4, SCALAR);
  set_variable_equation_type(4, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(4, "cM, grad(cM), psi, grad(psi), irxn");
  set_dependencies_gradient_term_RHS(4, "cM, grad(cM), grad(Phi)");

  // Variable 5
  set_variable_name(5, "cP");
  set_variable_type(5, SCALAR);
  set_variable_equation_type(5, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(5, "cP, grad(cP), psi, grad(psi)");
  set_dependencies_gradient_term_RHS(5, "cP, grad(cP), grad(Phi)");

  // Variable 6
  set_variable_name(6, "Phi");
  set_variable_type(6, SCALAR);
  set_variable_equation_type(6, TIME_INDEPENDENT);

  set_dependencies_value_term_LHS(
    6,
    "n, psi, grad(psi), cM, grad(cM), cP, Phi, grad(Phi), change(Phi), irxn");
  set_dependencies_gradient_term_LHS(6, "n, psi, cM, cP, grad(change(Phi))");
  set_dependencies_value_term_RHS(6, "grad(psi), irxn");
  set_dependencies_gradient_term_RHS(6, "psi, grad(Phi), grad(cM), grad(cP)");

  // Variable 7
  set_variable_name(7, "irxn");
  set_variable_type(7, SCALAR);
  set_variable_equation_type(7, AUXILIARY);

  set_dependencies_value_term_RHS(7, "cM, cP, Phi");
  set_dependencies_gradient_term_RHS(7, "");
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

  // The order parameter and its derivatives
  scalarvalueType n = variable_list.get_scalar_value(0);

  // The chemical potential and its derivatives
  scalargradType mux = variable_list.get_scalar_gradient(1);

  // The domain parameter and its derivatives
  scalarvalueType psi  = variable_list.get_scalar_value(2);
  scalargradType  psix = variable_list.get_scalar_gradient(2);

  // The chemical potential of the domain patameter
  scalargradType mupsix = variable_list.get_scalar_gradient(3);

  // The concentration of metal ion and its derivatives
  scalarvalueType cM  = variable_list.get_scalar_value(4);
  scalargradType  cMx = variable_list.get_scalar_gradient(4);

  // The concentration of supporting cation and its derivatives
  scalarvalueType cP  = variable_list.get_scalar_value(5);
  scalargradType  cPx = variable_list.get_scalar_gradient(5);

  // The electrolite potential and its derivatives
  scalargradType Phix = variable_list.get_scalar_gradient(6);

  // The reaction current
  scalarvalueType irxn = variable_list.get_scalar_value(7);

  // Calculation of capped fields
  psi = std::min(psi, constV(1.0));
  psi = std::max(psi, constV(lthresh));
  n   = std::min(n, constV(1.0));
  n   = std::max(n, constV(lthresh));

  // --- Calculation of terms needed in multiple expressions  ---
  // Magnifude of the gradient of psi
  scalarvalueType magpsix = constV(0.0);
  for (int i = 0; i < dim; i++)
    {
      magpsix = magpsix + psix[i] * psix[i];
    }
  magpsix = std::sqrt(magpsix);
  // Inverse of psi
  scalarvalueType invpsi = constV(1.0) / psi;
  // Velocity of the interface (scalar)
  scalarvalueType v = -constV(VMV / (zMV * FarC)) * irxn;
  // The mobility (including dependence of psi
  scalarvalueType MnV = MconstV * psi * std::abs(irxn);
  // Products needed to calculate the concentraion residual terms
  // Miscellaneous dot products
  scalarvalueType psixcMx    = constV(0.0);
  scalarvalueType psixcPx    = constV(0.0);
  scalarvalueType psixcMPhix = constV(0.0);
  scalarvalueType psixcPPhix = constV(0.0);

  for (int i = 0; i < dim; i++)
    {
      psixcMx    = psixcMx + psix[i] * cMx[i];
      psixcPx    = psixcPx + psix[i] * cPx[i];
      psixcMPhix = psixcMPhix + psix[i] * cM * Phix[i];
      psixcPPhix = psixcPPhix + psix[i] * cP * Phix[i];
    }

  // --- Calculation of residual terms for n  ---
  scalarvalueType rnV  = n + v * delt * magpsix;
  scalargradType  rnxV = -MnV * delt * mux;

  // --- Calculation of residual terms for psi ---
  scalarvalueType rpsiV  = psi - v * delt * magpsix;
  scalargradType  rpsixV = -MnV * delt * mupsix;

  // --- Calculation of residual terms for cM ---
  scalarvalueType rcMV =
    cM + delt * (constV(DMV) * invpsi * psixcMx +
                 constV(DMV * zMV * FarC / (RV * TV)) * psixcMPhix * invpsi +
                 constV(1.0 / (zMV * FarC)) * invpsi * magpsix * irxn);
  scalargradType rcMxV =
    -delt * (DMV * cMx + constV(DMV * zMV * FarC / (RV * TV)) * cM * Phix);

  // --- Calculation of residual terms for cP ---
  scalarvalueType rcPV =
    cP + delt * (constV(DPV) * invpsi * psixcPx +
                 constV(DPV * zPV * FarC / (RV * TV)) * psixcPPhix * invpsi);
  scalargradType rcPxV =
    -delt * (DPV * cPx + constV(DPV * zPV * FarC / (RV * TV)) * cP * Phix);

  // --- Submitting the terms for the governing equations ---
  // Residuals terms for the equation to evolve the order parameter
  variable_list.set_scalar_value_term_RHS(0, rnV);
  variable_list.set_scalar_gradient_term_RHS(0, rnxV);

  // Residuals terms for the equation to evolve the domain parameter
  variable_list.set_scalar_value_term_RHS(2, rpsiV);
  variable_list.set_scalar_gradient_term_RHS(2, rpsixV);

  // Residuals for the equation to evolve the concentration of metal ion
  variable_list.set_scalar_value_term_RHS(4, rcMV);
  variable_list.set_scalar_gradient_term_RHS(4, rcMxV);

  // Residuals for the equation to evolve the concentration of supporting cation
  variable_list.set_scalar_value_term_RHS(5, rcPV);
  variable_list.set_scalar_gradient_term_RHS(5, rcPxV);
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

  // The order parameter and its derivatives
  scalarvalueType n  = variable_list.get_scalar_value(0);
  scalargradType  nx = variable_list.get_scalar_gradient(0);

  // The domain parameter and its derivatives
  scalarvalueType psi  = variable_list.get_scalar_value(2);
  scalargradType  psix = variable_list.get_scalar_gradient(2);

  // The concentration of metal ion and its derivatives
  scalarvalueType cM  = variable_list.get_scalar_value(4);
  scalargradType  cMx = variable_list.get_scalar_gradient(4);

  // The concentration of supporting cation and its derivatives
  scalarvalueType cP  = variable_list.get_scalar_value(5);
  scalargradType  cPx = variable_list.get_scalar_gradient(5);

  // The electrolite potential and its derivatives
  scalarvalueType Phi  = variable_list.get_scalar_value(6);
  scalargradType  Phix = variable_list.get_scalar_gradient(6);

  // The reaction current
  scalarvalueType irxn = variable_list.get_scalar_value(7);

  // Calculation of capped fields
  psi = std::min(psi, constV(1.0));
  psi = std::max(psi, constV(lthresh));
  n   = std::min(n, constV(1.0));
  n   = std::max(n, constV(lthresh));

  // --- Calculation of terms needed in multiple expressions  ---
  // Magnitude of the gradient of the domain parameter
  scalarvalueType magpsix = constV(0.0);
  for (int i = 0; i < dim; i++)
    {
      magpsix = magpsix + psix[i] * psix[i];
    }
  magpsix = std::sqrt(magpsix);
  // Normal vector
  scalargradType nvec = psix / (magpsix + lthresh);

  // Calculation of gradPhi dot n
  scalarvalueType gradPhin = constV(0.0);
  for (int i = 0; i < dim; i++)
    {
      gradPhin = gradPhin + Phix[i] * nvec[i];
    }

  // --- Calculation of residual terms for mu  ---
  // Derivative of bulk free energy with respect to n
  scalarvalueType fnV = WV * n * (n * n - constV(1.0) + constV(2.0) * gammaV * psi * psi);
  scalarvalueType rmuV  = fnV;
  scalargradType  rmuxV = epssqV * nx;

  // --- Calculation of residual terms for mupsi  ---
  // Derivative of bulk free energy with respect to psi
  scalarvalueType fpsiV =
    WV * psi * (psi * psi - constV(1.0) + constV(2.0) * gammaV * n * n);
  scalarvalueType rmupsiV  = fpsiV;
  scalargradType  rmupsixV = epssqV * psix;

  // --- Calculation of residual terms for irxn  ---
  // Overpotential
  scalarvalueType eta = VsV - EcorrV - Phi;
  // Maximum current prefactor
  scalarvalueType prefac = zMV * FarC / (constV(1.0) - VMV * cM);
  // 2*deltaV/tau
  scalarvalueType twodelintau = (zMV * DMV * FarC / (RV * TV)) * std::abs(gradPhin);
  for (unsigned int j = 0; j < psi.size(); j++)
    {
      if (twodelintau[j] < DMV / (2.0 * deltaV))
        twodelintau[j] = DMV / (2.0 * deltaV);
    }
  scalarvalueType term_1 = twodelintau * (constV(cMsatV) - cM);
  scalarvalueType term_2 = constV(0.0);
  for (int i = 0; i < dim; i++)
    {
      term_2 = term_2 + cMx[i] * nvec[i];
    }
  term_2                 = DMV * term_2;
  scalarvalueType term_3 = (zMV * FarC * DMV * cM / (RV * TV)) * gradPhin;
  scalarvalueType imax   = prefac * (term_1 + term_2 + term_3);
  // Capping imax
  imax = std::min(imax, constV(imax_max));
  imax = std::max(imax, constV(imax_min));

  // Exponential term
  scalarvalueType eterm = std::exp(constV(expconstV) * eta);

  // reaction current residual
  scalarvalueType rirxnV =
    constV(icorrV) * eterm / (constV(1.0) + constV(icorrV) * eterm / imax);

  // --- Calculation of residual terms for Phi  ---
  scalarvalueType kappa =
    constV(FarC * FarC / (RV * TV)) * (constV(zMV * (zMV * DMV - znV * DnV)) * cM +
                                       constV(zPV * (zPV * DPV - znV * DnV)) * cP);
  scalarvalueType Msum = constV(FarC * zMV * (DnV - DMV));
  scalarvalueType Psum = constV(FarC * zPV * (DnV - DPV));

  scalarvalueType rPhiV  = -magpsix * irxn;
  scalargradType  rPhixV = psi * kappa * Phix - Msum * psi * cMx - Psum * psi * cPx;

  // --- Submitting the terms for the governing equations ---
  // Residuals for the equation to calculate mu
  variable_list.set_scalar_value_term_RHS(1, rmuV);
  variable_list.set_scalar_gradient_term_RHS(1, rmuxV);

  // Residuals for the equation to calculate mupsi
  variable_list.set_scalar_value_term_RHS(3, rmupsiV);
  variable_list.set_scalar_gradient_term_RHS(3, rmupsixV);

  // Residuals for the equation to evolve the Potential
  variable_list.set_scalar_value_term_RHS(6, rPhiV);
  variable_list.set_scalar_gradient_term_RHS(6, rPhixV);

  // Residuals for the equation to evolve irxn
  variable_list.set_scalar_value_term_RHS(7, rirxnV);
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
  // The order parameter and its derivatives
  scalarvalueType n = variable_list.get_scalar_value(0);

  // The domain parameter and its derivatives
  scalarvalueType psi  = variable_list.get_scalar_value(2);
  scalargradType  psix = variable_list.get_scalar_gradient(2);

  // The concentration of metal ion and its derivatives
  scalarvalueType cM  = variable_list.get_scalar_value(4);
  scalargradType  cMx = variable_list.get_scalar_gradient(4);

  // The concentration of supporting cation and its derivatives
  scalarvalueType cP = variable_list.get_scalar_value(5);

  // The electrolite potential and its derivatives
  scalarvalueType Phi  = variable_list.get_scalar_value(6);
  scalargradType  Phix = variable_list.get_scalar_gradient(6);

  // The change in potential in the electrode and its derivatives
  scalarvalueType DPhi  = variable_list.get_change_in_scalar_value(6);
  scalargradType  DPhix = variable_list.get_change_in_scalar_gradient(6);

  // Calculation of capped fields
  psi = std::min(psi, constV(1.0));
  psi = std::max(psi, constV(lthresh));
  n   = std::min(n, constV(1.0));
  n   = std::max(n, constV(lthresh));

  // --- Calculation of terms needed in multiple expressions  ---
  // Magnitude of the gradient of the domain parameter
  scalarvalueType magpsix = constV(0.0);
  for (int i = 0; i < dim; i++)
    {
      magpsix = magpsix + psix[i] * psix[i];
    }
  magpsix = std::sqrt(magpsix);

  // Normal vector
  scalargradType nvec = psix / (magpsix + lthresh);

  // Calculation of gradPhi dot n
  scalarvalueType gradPhin = constV(0.0);
  for (int i = 0; i < dim; i++)
    {
      gradPhin = gradPhin + Phix[i] * nvec[i];
    }

  // --- Calculation of residual terms for DPhi  ---
  scalarvalueType kappa =
    constV(FarC * FarC / (RV * TV)) * (constV(zMV * (zMV * DMV - znV * DnV)) * cM +
                                       constV(zPV * (zPV * DPV - znV * DnV)) * cP);

  // --- Calculation of residual terms for irxnp (the derivative of irxn wrt to
  // Phi)   --- Overpotential
  scalarvalueType eta = VsV - EcorrV - Phi;
  // Maximum current prefactor
  scalarvalueType prefac = zMV * FarC / (constV(1.0) - VMV * cM);
  // 2*deltaV/tau

  scalarvalueType twodelintau = (zMV * DMV * FarC / (RV * TV)) * std::abs(gradPhin);
  for (unsigned int j = 0; j < psi.size(); j++)
    {
      if (twodelintau[j] < DMV / (2.0 * deltaV))
        twodelintau[j] = DMV / (2.0 * deltaV);
    }
  scalarvalueType term_1 = twodelintau * (cMsatV - cM);
  scalarvalueType term_2 = constV(0.0);
  for (int i = 0; i < dim; i++)
    {
      term_2 = term_2 + cMx[i] * nvec[i];
    }
  term_2                 = DMV * term_2;
  scalarvalueType term_3 = (zMV * FarC * DMV * cM / (RV * TV)) * gradPhin;
  scalarvalueType imax   = prefac * (term_1 + term_2 + term_3);
  // Capping imax
  imax = std::min(imax, constV(imax_max));
  imax = std::max(imax, constV(imax_min));

  // Exponential term
  scalarvalueType eterm   = std::exp(constV(expconstV) * eta);
  scalarvalueType numer_i = icorrV * eterm;
  scalarvalueType denom_i = constV(1.0) + icorrV * eterm / imax;
  scalarvalueType xip     = -constV(expconstV);
  // scalarvalueType irxnp = (denom_i*icorrV*xip*eterm -
  // numer_i*icorrV*xip*eterm/imax)/(denom_i*denom_i);
  scalarvalueType irxnp =
    icorrV * xip * eterm * (denom_i - numer_i / imax) / (denom_i * denom_i);

  scalarvalueType rDPhi  = magpsix * irxnp * DPhi;
  scalargradType  rDPhix = -psi * kappa * DPhix;

  // Residuals for the equation to evolve the order parameter
  variable_list.set_scalar_value_term_LHS(6, rDPhi);
  variable_list.set_scalar_gradient_term_LHS(6, rDPhix);
}

// =================================================================================
// thresholdField: a function particular to this app
// =================================================================================
// Method that caps the value of the order parameter and the domain parameter
template <int dim, int degree>
void
customPDE<dim, degree>::capFields(VectorizedArray<double> &ncp,
                                  VectorizedArray<double> &psicp,
                                  VectorizedArray<double>  n,
                                  VectorizedArray<double>  psi) const
{
  // Capping n to lower threshold bound and upper bound of 1
  for (unsigned j = 0; j < ncp.size(); j++)
    {
      ncp[j] = n[j];
      if (n[j] < 0.0)
        ncp[j] = 0.0;
      if (n[j] > 1.0)
        ncp[j] = 1.0;
    }
  // Capping psi to lower threshold bound and upper bound of 1
  for (unsigned j = 0; j < ncp.size(); j++)
    {
      psicp[j] = psi[j];
      if (psi[j] < lthresh)
        psicp[j] = lthresh;
      if (psi[j] > 1.0)
        psicp[j] = 1.0;
    }
}
