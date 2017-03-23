//initial condition
template <int dim>
double InitialCondition<dim>::value (const dealii::Point<dim> &p, const unsigned int component) const
{
  double scalar_IC=0;
	  // =====================================================================
	  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
	  // =====================================================================
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  // Initial condition parameters
	  #define x_denom (5.0/2.0*scaleFactor)*(5.0/2.0*scaleFactor)
	  #define y_denom (5.0/2.0*scaleFactor)*(5.0/2.0*scaleFactor)
	  #define z_denom (20.0/2.0*scaleFactor)*(20.0/2.0*scaleFactor)
	  #define initial_interface_coeff (0.01*scaleFactor)
	  #define initial_radius 1.0
	  #define c_matrix 1.0e-6
      #define c_precip 0.125
	  #define adjust_avg_c false
	  #define c_avg 0.004

	  //set result equal to the structural order parameter initial condition
	  double dx=spanX/( (double)subdivisionsX )/std::pow(2.0,refineFactor);
	  double dy=spanY/( (double)subdivisionsY )/std::pow(2.0,refineFactor);
	  double dz=spanZ/( (double)subdivisionsZ )/std::pow(2.0,refineFactor);
	  double r=0.0;
	  std::vector<double> ellipsoid_denoms;
	  ellipsoid_denoms.push_back(x_denom);
	  ellipsoid_denoms.push_back(y_denom);
	  ellipsoid_denoms.push_back(z_denom);

	  double precip_center[3] = {0.0*scaleFactor, 0.0*scaleFactor, 0.0*scaleFactor};

	  for (unsigned int i=0; i<dim; i++){
		  r += (p.operator()(i)-precip_center[i])*(p.operator()(i)-precip_center[i])/ellipsoid_denoms[i];
	  }
	  r = sqrt(r);


	  if (index==0){
		  scalar_IC = 0.5*(c_precip-c_matrix)*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff))) + c_matrix;
	  }
	  else if (index==1){
		  scalar_IC = 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));
	  }

	  // =====================================================================
	  return scalar_IC;
}

template <int dim, int degree>
void customPDE<dim,degree>::shiftConcentration(unsigned int index){
#ifdef adjust_avg_c

	double integrated_concentration;
	this->computeIntegral(integrated_concentration,index);

	double volume = spanX;
	if (dim > 1) {
		volume *= spanY;
		if (dim > 2) {
			volume *= spanZ;
		}
	}

	double shift = c_avg - integrated_concentration/volume;

	if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
		std::cout<<"Matrix concentration shifted from " <<c_matrix<<" to " << c_matrix+shift <<std::endl;
	}

	try{
		if (shift + c_matrix < 0.0) {throw 0;}
	}
	catch (int e){
		Assert (shift > c_matrix, ExcMessage("An exception occurred. Initial concentration was shifted below zero."));
	}

	for (unsigned int dof=0; dof<this->solutionSet[index]->local_size(); ++dof){
		this->solutionSet[index]->local_element(dof) += shift;
	}

	this->solutionSet[index]->update_ghost_values();

	this->computeIntegral(integrated_concentration,index);

#endif
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

	  if (index==2){
		  vector_IC(0) = 0.0;
		  vector_IC(1) = 0.0;
          if (dim == 3){
        	  vector_IC(2) = 0.0;
          }
	  }
	  // =====================================================================
}

template <int dim, int degree>
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
	// Third input: BC type (options are "ZERO_DERIVATIVE" and "DIRICHLET")
	// Fourth input: BC value (ignored unless the BC type is "DIRICHLET")
	// Odd inputs after the third: BC type
	// Even inputs after the third: BC value
	// Face numbering: starts at zero with the minimum of the first direction, one for the maximum of the first direction
	//						two for the minimum of the second direction, etc.

	this->inputBCs(0,0,"ZERO_DERIVATIVE",0);

	this->inputBCs(1,0,"ZERO_DERIVATIVE",0);

	if (dim == 2){
		this->inputBCs(2,0,"DIRICHLET",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0);
		this->inputBCs(2,1,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0);
	}
	else if (dim == 3){
		this->inputBCs(2,0,"DIRICHLET",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0);
		this->inputBCs(2,1,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0);
		this->inputBCs(2,2,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0);
	}

}


