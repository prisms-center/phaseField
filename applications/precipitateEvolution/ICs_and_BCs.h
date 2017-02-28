//initial condition
template <int dim>
double InitialCondition<dim>::value (const dealii::Point<dim> &p, const unsigned int component) const
 {
  double scalar_IC=0.0;
	  // =====================================================================
	  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
	  // =====================================================================
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  double dx=spanX/((double) subdivisionsX)/std::pow(2.0,refineFactor);
	  double r=0.0;

	  if (index==0){
		  scalar_IC = 0.04;
	  }
	  else if (index==1){
		  #if problemDIM == 2
		  double r1=p.distance(Point<dim>(spanX/3.0,spanY/3.0));
		  double r2=p.distance(Point<dim>(2*spanX/3.0,2*spanY/3.0));
		  #elif problemDIM == 3
		  double r1=p.distance(Point<dim>(spanX/3.0,spanY/3.0,0.5*spanZ));
		  double r2=p.distance(Point<dim>(2*spanX/3.0,2*spanY/3.0,0.5*spanZ));
		  #endif
		  r=std::min(r1,r2);
		  scalar_IC = 0.5*(1.0-std::tanh((r-spanX/16.0)/(dx)));
	  }
	  else if (index==2){
		  #if problemDIM == 2
		  r=p.distance(Point<dim>(3*spanX/4.0,spanY/4.0));
		  #elif problemDIM == 3
		  r=p.distance(Point<dim>(3*spanX/4.0,spanY/4.0,0.5*spanZ));
		  #endif
	        scalar_IC = 0.5*(1.0-std::tanh((r-spanX/16.0)/(dx)));
	  }
	  else if (index==3){
		  #if problemDIM == 2
		  r=p.distance(Point<dim>(spanX/4.0,3.0*spanY/4.0));
		  #elif problemDIM == 3
		  r=p.distance(Point<dim>(spanX/4.0,3.0*spanY/4.0,0.5*spanZ));
		  #endif
		  scalar_IC = 0.5*(1.0-std::tanh((r-spanX/16.0)/(dx)));
	  }

	  // =====================================================================
	  return scalar_IC;
  }

template <int dim>
 void InitialConditionVec<dim>::vector_value (const dealii::Point<dim> &p, dealii::Vector<double> &vector_IC) const
 {
	  // =====================================================================
	  // ENTER THE INITIAL CONDITIONS HERE FOR VECTOR FIELDS
	  // =====================================================================
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  if (index==4){
		  vector_IC(0) = 0.0;
		  vector_IC(1) = 0.0;
          if (dim == 3){
        	  vector_IC(2) = 0.0;
          }
	  }
	  // =====================================================================
  }

 template class InitialCondition<2>;
 template class InitialConditionVec<2>;
 template class InitialCondition<3>;
 template class InitialConditionVec<3>;

template <int dim,int degree>
void customPDE<dim,degree>::setBCs(){

	// =====================================================================
	// ENTER THE BOUNDARY CONDITIONS HERE
	// =====================================================================
	// This function sets the BCs for the problem variables
	// The function "inputBCs" should be called for each component of
	// each variable and should be in numerical order. Four input arguments
	// set the same BC on the entire boundary. Two plus two times the
	// number of dimensions inputs sets separate BCs on each face of the domain.
	// Inputs to "inputBCs":
	// First input: variable number
	// Second input: component number
	// Third input: BC type (options are "ZERO_DERIVATIVE", "DIRICHLET", and "PERIODIC")
	// Fourth input: BC value (ignored unless the BC type is "DIRICHLET")
	// Odd inputs after the third: BC type
	// Even inputs after the third: BC value
	// Face numbering: starts at zero with the minimum of the first direction, one for the maximum of the first direction
	//						two for the minimum of the second direction, etc.

	this->inputBCs(0,0,"ZERO_DERIVATIVE",0);

	this->inputBCs(1,0,"ZERO_DERIVATIVE",0);

	this->inputBCs(2,0,"ZERO_DERIVATIVE",0);

	this->inputBCs(3,0,"ZERO_DERIVATIVE",0);

	this->inputBCs(4,0,"DIRICHLET",0.0);
	this->inputBCs(4,1,"DIRICHLET",0.0);
	if (dim == 3){
		this->inputBCs(4,2,"DIRICHLET",0.0);
	}

}


