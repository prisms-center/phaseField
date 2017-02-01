// List of variables and residual equations for the coupled Allen-Cahn/Cahn-Hilliard with regularized anisotropic interfacial energy example application

// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 3

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"c", "n","biharm"}
#define variable_type {"SCALAR","SCALAR","SCALAR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC","PARABOLIC"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define need_val {true, true, true}
#define need_grad {true, true, true}
#define need_hess {false, false, false}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_val_residual {true, true, false}
#define need_grad_residual {true, true, true}

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

//define Cahn-Hilliard parameters (No Gradient energy term)
#define McV 1.0

//define free energy expressions
#define faV (0.5*c*c/16.0)
#define facV (0.5*c/8.0)
#define faccV (0.5/8.0)
#define fbV (0.5*(c-1.0)*(c-1.0)/16.0)
#define fbcV (0.5*(c-1.0)/8.0)
#define fbccV (0.5/8.0)
#define hV (3.0*n*n-2.0*n*n*n)
#define hnV (6.0*n-6.0*n*n)

//anisotropy and regularization parameters
#define epsilonM 0.2
#define delta2 1.0

//anisotropy gamma as a function of the components of the normal vector
//current anisotropy has 4-fold or octahedral symmetry
#if problemDIM==1
#define gamma 1.0
#elif problemDIM==2
//writing out powers instead of using std::pow(double,double) for performance reasons
#define gamma (1.0+epsilonM*(4.0*(normal[0]*normal[0]*normal[0]*normal[0]+normal[1]*normal[1]*normal[1]*normal[1])-3.0))
#else
#define gamma (1.0+epsilonM*(4.0*(normal[0]*normal[0]*normal[0]*normal[0]+normal[1]*normal[1]*normal[1]*normal[1]+normal[2]*normal[2]*normal[2]*normal[2])-3.0))
#endif

//derivatives of gamma with respect to the components of the unit normal
#define gammanx (epsilonM*16.0*normal[0]*normal[0]*normal[0])
#define gammany (epsilonM*16.0*normal[1]*normal[1]*normal[1])
#define gammanz (epsilonM*16.0*normal[2]*normal[2]*normal[2])

//Allen-Cahn mobility (isotropic)
#define MnV 0.1

//Allen-Cahn mobility (anisotropic)
//#define MnV (1.0/(gamma*gamma+1e-10))

//define required residuals (aniso defined in model)
#define rbiharmV constV(0.0)
#define rbiharmxV (-nx)
#define rcV   ( c )
#define rcxV  (constV(-McV*timeStep)*(cx*((1.0-hV)*faccV+hV*fbccV)+nx*hnV*(fbcV-facV)))
#define rnV  (n-constV(timeStep*MnV)*(fbV-faV)*hnV)
#define rnxV (constV(timeStep*MnV)*(-aniso+constV(delta2)*biharmx))

// =================================================================================
// residualRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "modelVariablesList" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of
// that quadrature point is given by "q_point_loc". The function outputs
// "modelResidualsList", a list of the value and gradient terms of the residual for
// each residual equation. The index for each variable in these lists corresponds to
// the order it is defined at the top of this file (starting at 0).template <int dim>
template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<modelVariable<dim> > & modelVariablesList,
     std::vector<modelResidual<dim> > & modelResidualsList,
     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//c
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

//n
scalarvalueType n = modelVariablesList[1].scalarValue;
scalargradType nx = modelVariablesList[1].scalarGrad;

//biharm
scalargradType biharmx = modelVariablesList[2].scalarGrad;

// anisotropy code
scalarvalueType normgradn = std::sqrt(nx.norm_square());
scalargradType normal = nx/(normgradn+constV(1.0e-16));
scalarvalueType gamma_scl = gamma;
scalargradType aniso;
#if problemDIM==1
      aniso = gamma_scl*gamma_scl*nx
#else
      scalargradType dgammadnorm;
      dgammadnorm[0]=gammanx;
      dgammadnorm[1]=gammany;
#if problemDIM>2
      dgammadnorm[2]=gammanz;
#endif
      for (unsigned int i=0; i<problemDIM; ++i){
	      for (unsigned int j=0; j<problemDIM; ++j){
		      aniso[i] += -normal[i]*normal[j]*dgammadnorm[j];
		      if (i==j) aniso[i] +=dgammadnorm[j];
	      }
      }
      aniso = gamma_scl*(aniso*normgradn+gamma_scl*nx);
#endif
// end anisotropy code

modelResidualsList[0].scalarValueResidual = rcV;
modelResidualsList[0].scalarGradResidual = rcxV;

modelResidualsList[1].scalarValueResidual = rnV;
modelResidualsList[1].scalarGradResidual = rnxV;

modelResidualsList[2].scalarGradResidual = rbiharmxV;
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
template <int dim>
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
    modelResidual<dim> & modelRes, dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

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
template <int dim>
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim> > & modelVarList,
    const dealii::VectorizedArray<double> & JxW_value,
    dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {

scalarvalueType total_energy_density = constV(0.0);

//c
scalarvalueType c = modelVarList[0].scalarValue;
scalargradType cx = modelVarList[0].scalarGrad;

//n1
scalarvalueType n = modelVarList[1].scalarValue;
scalargradType nx = modelVarList[1].scalarGrad;

//biharm
scalarvalueType biharm = modelVarList[2].scalarValue;

scalarvalueType f_chem = (constV(1.0)-hV)*faV + hV*fbV;

scalargradType normal = nx/(std::sqrt(nx.norm_square())+constV(1.0e-16));
scalarvalueType f_grad = constV(0.5)*gamma*gamma*nx*nx;
scalarvalueType f_reg = constV(0.5*delta2)*biharm*biharm;

total_energy_density = f_chem + f_grad + f_reg;

assembler_lock.acquire ();
for (unsigned i=0; i<c.n_array_elements;i++){
  if (c[i] > 1.0e-10){
	  this->energy+=total_energy_density[i]*JxW_value[i];
	  this->energy_components[0]+= f_chem[i]*JxW_value[i];
	  this->energy_components[1]+= f_grad[i]*JxW_value[i];
	  this->energy_components[2]+= f_reg[i]*JxW_value[i];
  }
}
assembler_lock.release ();
}
