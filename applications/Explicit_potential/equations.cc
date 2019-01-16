// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"n");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "n, grad(psi), Phi");
    set_dependencies_gradient_term_RHS(0, "grad(mu)");

	// Variable 1
	set_variable_name				(1,"mu");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(1, "n, psi");
    set_dependencies_gradient_term_RHS(1, "grad(n)");

    // Variable 2
    set_variable_name				(2,"psi");
    set_variable_type				(2,SCALAR);
    set_variable_equation_type		(2,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(2, "n");

    // Variable 3
    set_variable_name				(3,"cM");
    set_variable_type				(3,SCALAR);
    set_variable_equation_type		(3,EXPLICIT_TIME_DEPENDENT);
    
    set_dependencies_value_term_RHS(3, "cM, grad(cM), psi, grad(psi)");
    set_dependencies_gradient_term_RHS(3, "cM, grad(cM), grad(Phi)");
    
    // Variable 4
    set_variable_name				(4,"cP");
    set_variable_type				(4,SCALAR);
    set_variable_equation_type		(4,EXPLICIT_TIME_DEPENDENT);
    
    set_dependencies_value_term_RHS(4, "cP, grad(cP), psi, grad(psi)");
    set_dependencies_gradient_term_RHS(4, "cP, grad(cP), grad(Phi)");
   
    // Variable 5 (Making Phi temporarily time-dependent)
    set_variable_name                (5,"Phi");
    set_variable_type                (5,SCALAR);
    set_variable_equation_type        (5,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(5, "Phi");

    // Variable 5
    /*
    set_variable_name				(5,"Phi");
    set_variable_type				(5,SCALAR);
    set_variable_equation_type		(5,TIME_INDEPENDENT);
    
    set_dependencies_value_term_RHS(5, "grad(psi), Phi");
    set_dependencies_gradient_term_RHS(5, "n, psi, grad(Phi), cM, cP, grad(cM), grad(cP)");
    set_dependencies_value_term_LHS(5, "grad(psi), cM, Phi, change(Phi)");
    set_dependencies_gradient_term_LHS(5, "n, psi, cM, cP, grad(change(Phi))");
	*/

}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

    //Timestep
    scalarvalueType delt = constV(userInputs.dtValue);
    
    // The order parameter and its derivatives
    scalarvalueType n = variable_list.get_scalar_value(0);
    scalargradType nx = variable_list.get_scalar_gradient(0);
    
    // The chemical potential and its derivatives
    scalargradType mux = variable_list.get_scalar_gradient(1);
    
    // The domain parameter and its derivatives
    scalarvalueType psi = variable_list.get_scalar_value(2);
    scalargradType psix = variable_list.get_scalar_gradient(2);
    
    // The concentration of metal ion and its derivatives
    scalarvalueType cM = variable_list.get_scalar_value(3);
    scalargradType cMx = variable_list.get_scalar_gradient(3);
    
    // The concentration of supporting cation and its derivatives
    scalarvalueType cP = variable_list.get_scalar_value(4);
    scalargradType cPx = variable_list.get_scalar_gradient(4);
    
    // The electrolite potential and its derivatives
    scalarvalueType Phi = variable_list.get_scalar_value(5);
    scalargradType Phix = variable_list.get_scalar_gradient(5);
    
    // Parameters in the explicit equations can be set here
    
    //Calculation of capped fields
    dealii::VectorizedArray<double> ncp = constV(0.0);
    dealii::VectorizedArray<double> psicp = constV(0.0);
    psi=1.0-n;
    capFields(ncp, psicp, n, psi);
    n = ncp;
    psi = psicp;
    
    // The derivative of the local free energy
    scalarvalueType fnV = n*(n*n-constV(1.0)+constV(3.0)*psi*psi);
    
    //Magnitude of the gradient of the domain parameter
    scalarvalueType magpsix = constV(0.0);
    for (int i=0; i<dim; i++){
        magpsix = magpsix + psix[i]*psix[i];
    }
    magpsix = std::sqrt(magpsix);
    
    //Miscellaneous dot products
    scalarvalueType psixcMx = constV(0.0);
    scalarvalueType psixcPx = constV(0.0);
    scalarvalueType psixcMPhix = constV(0.0);
    scalarvalueType psixcPPhix = constV(0.0);
    
    for (int i=0; i<dim; i++){
        psixcMx = psixcMx + psix[i]*cMx[i];
        psixcPx = psixcPx + psix[i]*cPx[i];
        psixcMPhix = psixcMPhix + psix[i]*cM*Phix[i];
        psixcPPhix = psixcPPhix + psix[i]*cP*Phix[i];
    }
    
    //Inverse of psi
    scalarvalueType invpsi=constV(1.0)/psi;
    
    //Overpotential
    scalarvalueType eta=VsV-EcorrV-Phi;
    
    //Reaction Current
    scalarvalueType irxn1=constV(i0V)*std::exp(constV(-zMV*betaV*FarC/(RV*TV))*eta);
    scalarvalueType irxn2=constV(i0V)*(constV(1.0)-cM/constV(cMsatV))*std::exp(constV(zMV*(1.0-betaV)*FarC/(RV*TV))*eta);
    scalarvalueType irxn = irxn1+irxn2;
    
    //Velocity Field
    scalarvalueType v=-constV(VMV/(zMV*FarC))*irxn;
    
    // The residuals
    scalarvalueType rnV = n + v*delt*magpsix;
    scalargradType rnxV = constV(-MnV)*psi*delt*mux;
    
    scalarvalueType rmuV = fnV;
    scalargradType rmuxV = constV(epssqV)*nx;
    
    scalarvalueType rpsiV = psi;
    scalargradType rpsixV = psix;
    
    scalarvalueType rcMV = cM + delt*(constV(DMV)*invpsi*psixcMx + constV(DMV*zMV*FarC/(RV*TV))*psixcMPhix*invpsi + constV(1.0/(zMV*FarC))*invpsi*magpsix*irxn);
    scalargradType rcMxV = -delt*(DMV*cMx + constV(DMV*zMV*FarC/(RV*TV))*cM*Phix);
    
    scalarvalueType rcPV = cP + delt*(constV(DPV)*invpsi*psixcPx + constV(DPV*zPV*FarC/(RV*TV))*psixcPPhix*invpsi);
    scalargradType rcPxV = -delt*(DPV*cPx + constV(DPV*zPV*FarC/(RV*TV))*cP*Phix);

    scalarvalueType rPhiV = Phi;
    //scalargradType rPhixV = constV(0.0)*Phix;

    
    // Residuals for the equation to evolve the order parameter
    variable_list.set_scalar_value_term_RHS(0,rnV);
    variable_list.set_scalar_gradient_term_RHS(0,rnxV);
    
    // Residuals for the equation to evolve the chemical potential
    variable_list.set_scalar_value_term_RHS(1,rmuV);
    variable_list.set_scalar_gradient_term_RHS(1,rmuxV);
    
    // Residuals for the equation to evolve the domain parameter
    variable_list.set_scalar_value_term_RHS(2,rpsiV);
    variable_list.set_scalar_gradient_term_RHS(2,rpsixV);
    
    // Residuals for the equation to evolve the concentration of metal ion
    variable_list.set_scalar_value_term_RHS(3,rcMV);
    variable_list.set_scalar_gradient_term_RHS(3,rcMxV);
    
    // Residuals for the equation to evolve the concentration of supporting cation
    variable_list.set_scalar_value_term_RHS(4,rcPV);
    variable_list.set_scalar_gradient_term_RHS(4,rcPxV);
    
    variable_list.set_scalar_value_term_RHS(5,rPhiV);
    //variable_list.set_scalar_gradient_term_RHS(5,rPhixV);

}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}

// =================================================================================
// thresholdField: a function particular to this app
// =================================================================================
//Method that caps the value of the order parameter and the domain parameter
template <int dim,int degree>
void customPDE<dim,degree>::capFields(dealii::VectorizedArray<double> & ncp, dealii::VectorizedArray<double> & psicp,
                                      dealii::VectorizedArray<double> n, dealii::VectorizedArray<double> psi) const {
    //Capping n to lower threshold bound and upper bound of 1
    for (unsigned j=0; j<ncp.n_array_elements;j++){
        ncp[j]=n[j];
        if (n[j] < lthresh)
            ncp[j] = lthresh;
        if (n[j] > 1.0)
            ncp[j] = 1.0;
    }
    //Capping psi to lower threshold bound and upper bound of 1
    for (unsigned j=0; j<ncp.n_array_elements;j++){
        psicp[j]=psi[j];
        if (psi[j] < lthresh)
            psicp[j] = lthresh;
        if (psi[j] > 1.0)
            psicp[j] = 1.0;
    }
}
