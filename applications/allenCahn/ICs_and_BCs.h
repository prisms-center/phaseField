template <int dim>
double InitialCondition<dim>::value (const dealii::Point<dim> &p, const unsigned int component) const
 {
  double scalar_IC;
  // =====================================================================
  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
  // =====================================================================
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the initial condition for each variable
  // according to its variable index.

  double x_loc[12] = {0.1, 0.8, 0.5, 0.4, 0.3, 0.8, 0.9, 0.0, 0.1, 0.5, 1, 0.7};
  double y_loc[12] = {0.3, 0.7, 0.2, 0.4, 0.9, 0.1, 0.5, 0.1, 0.6, 0.6, 1, 0.95};
  double rad[12] =   {12, 14, 19, 16, 11, 12, 17, 15, 20, 10, 11, 14};
  double dist;
  scalar_IC = 0;
  for (unsigned int i=0; i<12; i++){
	  #if problemDIM == 2
		  dist = p.distance(dealii::Point<dim>(x_loc[i]*spanX,y_loc[i]*spanY));
	  #elif problemDIM == 3
		  dist = p.distance(dealii::Point<dim>(x_loc[i]*spanX,y_loc[i]*spanY,0.5*spanZ));
	  #endif
	  if (dist < rad[i]){
		  scalar_IC = 1.0;
	  }

  };

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


	  // =====================================================================
 }

 template class InitialCondition<2>;
 template class InitialConditionVec<2>;
 template class InitialCondition<3>;
 template class InitialConditionVec<3>;



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
	// Third input: BC type (options are "ZERO_DERIVATIVE", "DIRICHLET", and "PERIODIC")
	// Fourth input: BC value (ignored unless the BC type is "DIRICHLET")
	// Odd inputs after the third: BC type
	// Even inputs after the third: BC value
	// Face numbering: starts at zero with the minimum of the first direction, one for the maximum of the first direction
	//						two for the minimum of the second direction, etc.

	this->inputBCs(0,0,"ZERO_DERIVATIVE",0);

	// =====================================================================

}


