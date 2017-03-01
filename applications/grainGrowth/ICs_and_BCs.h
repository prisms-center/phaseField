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

	  double dx=spanX/((double) subdivisionsX)/std::pow(2.0,refineFactor);
	  double x_loc[10] = {0.3, 0.7, 0.5, 0.4, 0.3, 0.1, 0.1, 0.6, 0.7, 0.7};
	  double y_loc[10] = {0.3, 0.7, 0.1, 0.5, 0.9, 0.1, 0.7, 0.6, 0.4, 0.2};
	  double rad[10] =   {6, 7, 9, 8, 5, 6, 4, 3, 8, 5};
	  double dist;
	  scalar_IC = 0;

	  #if problemDIM == 2
	  	  dist = p.distance(Point<dim>(x_loc[index]*spanX,y_loc[index]*spanY));
	  #elif problemDIM == 3
	  	  dist = p.distance(Point<dim>(x_loc[index]*spanX,y_loc[index]*spanY,0.5*spanZ));
	  #endif


	  scalar_IC = 0.5*(1.0-std::tanh((dist-rad[index])/(dx)));

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

	this->inputBCs(0,0,"PERIODIC",0);
	this->inputBCs(1,0,"PERIODIC",0);
	this->inputBCs(2,0,"PERIODIC",0);
	this->inputBCs(3,0,"PERIODIC",0);
	this->inputBCs(4,0,"PERIODIC",0);
	this->inputBCs(5,0,"PERIODIC",0);
	this->inputBCs(6,0,"PERIODIC",0);
	this->inputBCs(7,0,"PERIODIC",0);
	this->inputBCs(8,0,"PERIODIC",0);
	this->inputBCs(9,0,"PERIODIC",0);

	// =====================================================================

}


