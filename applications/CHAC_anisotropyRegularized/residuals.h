// Definition of the variables in the model
#define num_var 3
#define variable_name {"c", "n","biharm"}
#define variable_type {"SCALAR","SCALAR","SCALAR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC","PARABOLIC"}
#define need_val {true, true, true}
#define need_grad {true, true, true}
#define need_hess {false, false, false}
#define need_val_residual {true, true, false}
#define need_grad_residual {true, true, true}

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
#define gamma (1.0+epsilonM*(4.0*(std::pow(normal[0],4.0)+std::pow(normal[1],4.0))-3.0))
#else
#define gamma (1.0+epsilonM*(4.0*(std::pow(normal[0],4.0)+std::pow(normal[1],4.0)+std::pow(normal[2],4.0))-3.0))
#endif

//derivatives of gamma with respect to the components of the unit normal
#define gammanx (epsilonM*16.0*std::pow(normal[0],3.0))
#define gammany (epsilonM*16.0*std::pow(normal[1],3.0))
#define gammanz (epsilonM*16.0*std::pow(normal[2],3.0))

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

// ---------------------------------------------

template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<modelVariable<dim>> & modelVariablesList, 
     std::vector<modelResidual<dim>> & modelResidualsList, 
     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//c
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

//n
scalarvalueType n = modelVariablesList[1].scalarValue;
scalargradType nx = modelVariablesList[1].scalarGrad;

//biharm
scalarvalueType biharm = modelVariablesList[2].scalarValue;
scalargradType biharmx = modelVariablesList[2].scalarGrad;

// anisotropy code
scalarvalueType normgradn = std::sqrt(nx.norm_square());
scalargradType normal = nx/(normgradn+constV(1.0e-16));
scalargradType aniso = gamma*gamma*nx;
#if problemDIM>1
      scalargradType dgammadnorm;
      // projection matrix
      vectorgradType proj;
      for (unsigned int i=0; i<problemDIM; ++i){
	      for (unsigned int j=0; j<problemDIM; ++j){
		      proj[i][j] = -normal[i]*normal[j];
		      if (i==j) proj[i][j]+=constV(1.0);
	      }
      }
      dgammadnorm[0]=gammanx;
      dgammadnorm[1]=gammany;
#if problemDIM>2
      dgammadnorm[2]=gammanz;
#endif
      aniso += gamma*normgradn*(proj*dgammadnorm);
#endif
// end anisotropy code

modelResidualsList[0].scalarValueResidual = rcV;
modelResidualsList[0].scalarGradResidual = rcxV;

modelResidualsList[1].scalarValueResidual = rnV;
modelResidualsList[1].scalarGradResidual = rnxV;

modelResidualsList[2].scalarGradResidual = rbiharmxV;
}

template <int dim>
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim>> & modelVarList,
    modelResidual<dim> & modelRes, dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

}

template <int dim>
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim>> & modelVarList,
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
scalarvalueType f_reg = constV(0.5)*delta2*biharm*biharm;

total_energy_density = f_chem + f_grad + f_reg;

assembler_lock.acquire ();
for (unsigned i=0; i<c.n_array_elements;i++){
  // For some reason, some of the values in this loop
  if (c[i] > 1.0e-10){
	  this->energy+=total_energy_density[i]*JxW_value[i];
	  this->energy_components[0]+= f_chem[i]*JxW_value[i];
	  this->energy_components[1]+= f_grad[i]*JxW_value[i];
	  this->energy_components[2]+= f_reg[i]*JxW_value[i];
  }
}
assembler_lock.release ();
}
