
// List of variables and residual equations for the coupled Allen-Cahn/Cahn-Hilliard example application


// =================================================================================
// Define the Model and residual equations
// =================================================================================
// The model free energy equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

// =================================================================================
// Set the nucleation parameters
// =================================================================================

// Nucleation radius (order parameter)
#define semiaxis_a 5.0
#define semiaxis_b 5.0
#define semiaxis_c 5.0

// Hold time for order parameter
#define t_hold 20.0

// Small constant for sign function
#define epsil 1.0e-7

// Minimum distance between nuclei
#define minDistBetweenNuclei (4.0*semiaxis_a)
#define maxOrderParameterNucleation 0.01

// Number of time steps between nucleation attempts
#define skipNucleationSteps 30

// radius for order parameter hold
std::vector<double> opfreeze_semiaxes {1.5*semiaxis_a,1.5*semiaxis_b,1.5*semiaxis_c};

//Minimum distance from the edges of the system where nucleation can occur
#define borderreg (2.0*semiaxis_a)

// Constants k1 and k2 for nucleation rate in the bulk
#define k1 498.866
#define k2 4.14465

// Incubation time constant
#define tau 500.0

// =================================================================================

// Free energy for each phase and their first and second derivatives
#define faV (A0+A2*(c_alpha-calmin)*(c_alpha-calmin))
#define facV (2.0*A2*(c_alpha-calmin))
#define faccV (2.0*A2)
#define fbV (B0+B2*(c_beta-cbtmin)*(c_beta-cbtmin))
#define fbcV (2.0*B2*(c_beta-cbtmin))
#define fbccV (2.0*B2)

// Interpolation function and its derivative
#define hV (3.0*n*n - 2.0*n*n*n)
#define hnV (6.0*n - 6.0*n*n)

// KKS model c_alpha and c_beta as a function of c and h
#define c_alpha ((B2*(c-cbtmin*hV) + A2*calmin*hV)/(A2*hV+B2*(1.0-hV)))
#define c_beta ((A2*(c-calmin*(1.0-hV))+B2*cbtmin*(1.0-hV))/(A2*hV+B2*(1.0-hV)))

// Double-Well function (can be used to tune the interfacial energy)
#define fbarrierV (n*n - 2.0*n*n*n + n*n*n*n)
#define fbarriernV (2.0*n - 6.0*n*n + 4.0*n*n*n)

// Residual equations
// For concentration
#define term_muxV (cx + (c_alpha - c_beta)*hnV*nx)
#define rcV   (c)
#define rcxV  (constV(-McV*userInputs.dtValue)*term_muxV)
//For order parameter (gamma is a variable order parameter mobility factor)
#define rnV   (n-constV(userInputs.dtValue*MnV)*gamma*((fbV-faV)*hnV - (c_beta-c_alpha)*fbcV*hnV + W_barrier*fbarriernV))
#define rnxV  (constV(-userInputs.dtValue*KnV*MnV)*gamma*nx)

// =================================================================================
// residualRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "modelVariablesList" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of
// that quadrature point is given by "q_point_loc".
// This function also calculates the factor (gamma) that multiplies the order parameter mobility
// during the hold time after each nucleus has been seeded.
// The function outputs
// "modelResidualsList", a list of the value and gradient terms of the residual for
// each residual equation. The index for each variable in these lists corresponds to
// the order it is defined at the top of this file (starting at 0).
template <int dim, int degree>
void customPDE<dim,degree>::residualRHS(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

double time = this->currentTime;
double dx=userInputs.domain_size[0]/std::pow(2.0,userInputs.refine_factor);

//Calculation of mobility factor, gamma
dealii::VectorizedArray<double> gamma = constV(1.0);
dealii::VectorizedArray<double> x= q_point_loc[0];
dealii::VectorizedArray<double> y= q_point_loc[1];
dealii::VectorizedArray<double> z;

if (dim ==3)
    z= q_point_loc[2];


// The concentration and its derivatives (names here should match those in the macros above)
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

// The order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n = modelVariablesList[1].scalarValue;
scalargradType nx = modelVariablesList[1].scalarGrad;

dealii::VectorizedArray<double> nucleation_source_term = constV(0.0);

for (typename std::vector<nucleus<dim> >::const_iterator thisNucleus=nuclei.begin(); thisNucleus!=nuclei.end(); ++thisNucleus){

	// Calculate the weighted distance function to the order parameter freeze boundary (weighted_dist = 1.0 on that boundary)
	dealii::VectorizedArray<double> weighted_dist = constV(0.0);
	for (unsigned int i=0; i<dim; i++){
        dealii::VectorizedArray<double> temp = (thisNucleus->center(i) - q_point_loc(i));
        bool periodic_i = (userInputs.BC_list[1].var_BC_type[2*i]==PERIODIC);
        if (periodic_i){
            double domsize_i =userInputs.domain_size[i];
            for (unsigned j=0; j<n.n_array_elements;j++)
        		temp[j]=temp[j]-round(temp[j]/domsize_i)*domsize_i;
        }
        temp=temp/opfreeze_semiaxes[i];
		weighted_dist += temp*temp;
	}

	// Seed a nucleus if it was added to the list of nuclei this time step
	if (thisNucleus->seedingTimestep == this->currentIncrement){
    	for (unsigned i=0; i<n.n_array_elements;i++){
    		if (weighted_dist[i] <= 1.0){
    			// Find the weighted distance to the outer edge of the nucleus and use it to calculate the order parameter source term
    			double r = 0.0;
    			double avg_semiaxis = 0.0;
    			for (unsigned int j=0; j<dim; j++){
                    double temp = (thisNucleus->center(j) - q_point_loc(j)[i]);
                    bool periodic_j = (userInputs.BC_list[1].var_BC_type[2*j]==PERIODIC);
                    if (periodic_j){
                        double domsize_j =userInputs.domain_size[j];
                        temp=temp-round(temp/domsize_j)*domsize_j;
                    }
                    temp=temp/thisNucleus->semiaxes[j];
    				r += temp*temp;
    				avg_semiaxis += thisNucleus->semiaxes[j];
    			}
    			r = sqrt(r);
    			avg_semiaxis /= dim;
    			nucleation_source_term[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/interface_coeff));
    		}
    	}
    }
	dealii::Point<dim> r=thisNucleus->center;

    double nucendtime = thisNucleus->seededTime + thisNucleus->seedingTime;
    dealii::VectorizedArray<double> spacearg=(std::sqrt(weighted_dist)-constV(1.0))/constV(dx);
    dealii::VectorizedArray<double> timearg=constV(time-nucendtime)/constV(userInputs.dtValue);
    dealii::VectorizedArray<double> spacefactor=constV(0.5)-constV(0.5)*spacearg/(std::abs(spacearg)+epsil);
    dealii::VectorizedArray<double> timefactor=constV(0.5)-constV(0.5)*timearg/(std::abs(timearg)+epsil);
    dealii::VectorizedArray<double> localgamma= constV(1.0)-spacefactor*timefactor;
    gamma=gamma*localgamma;
}

// Residuals for the equation to evolve the concentration (names here should match those in the macros above)
modelResidualsList[0].scalarValueResidual = rcV;
modelResidualsList[0].scalarGradResidual = rcxV;

// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
modelResidualsList[1].scalarValueResidual = rnV + nucleation_source_term;
modelResidualsList[1].scalarGradResidual = rnxV;

}

// =================================================================================
// residualLHS (needed only if at least one equation is elliptic)
// =================================================================================
// This function calculates the residual equations for the iterative solver for
// elliptic equations.for each variable. It takes "modelVariablesList" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is given
// by "q_point_loc". The function outputs "modelRes", the value and gradient terms of
// for the left-hand-side of the residual equation for the iterative solver. The
// index for each variable in these lists corresponds to the order it is defined at
// the top of this file (starting at 0), not counting variables that have
// "need_val_LHS", "need_grad_LHS", and "need_hess_LHS" all set to "false". If there
// are multiple elliptic equations, conditional statements should be used to ensure
// that the correct residual is being submitted. The index of the field being solved
// can be accessed by "this->currentFieldIndex".
template <int dim, int degree>
void customPDE<dim,degree>::residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
		modelResidual<dim> & modelRes,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}

// =================================================================================
// energyDensity (needed only if calcEnergy == true)
// =================================================================================
// This function integrates the free energy density across the computational domain.
// It takes "modelVariablesList" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. It also
// takes the mapped quadrature weight, "JxW_value", as an input. The (x,y,z) location
// of the quadrature point is given by "q_point_loc". The weighted value of the
// energy density is added to "energy" variable and the components of the energy
// density are added to the "energy_components" variable (index 0: chemical energy,
// index 1: gradient energy, index 2: elastic energy).
template <int dim, int degree>
void customPDE<dim,degree>::energyDensity(const std::vector<modelVariable<dim> > & modelVariablesList,
											const dealii::VectorizedArray<double> & JxW_value,
											dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {

// The concentration and its derivatives (names here should match those in the macros above)
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

// The order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n = modelVariablesList[1].scalarValue;
scalargradType nx = modelVariablesList[1].scalarGrad;

// The homogenous free energy
scalarvalueType f_chem = (constV(1.0)-hV)*faV + hV*fbV;

// The gradient free energy
scalarvalueType f_grad = constV(0.5*KnV)*nx*nx;

// The total free energy
scalarvalueType total_energy_density;
total_energy_density = f_chem + f_grad;

// Loop to step through each element of the vectorized arrays. Working with deal.ii
// developers to see if there is a more elegant way to do this.
this->assembler_lock.acquire ();
for (unsigned i=0; i<c.n_array_elements;i++){
  if (c[i] > 1.0e-10){
	  this->energy+=total_energy_density[i]*JxW_value[i];
	  this->energy_components[0]+= f_chem[i]*JxW_value[i];
	  this->energy_components[1]+= f_grad[i]*JxW_value[i];
  }
}
this->assembler_lock.release ();
}
