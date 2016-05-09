// Initial geometry
#define x_denom 3.0625
#define y_denom 64.0
#define z_denom 42.25
#define initial_interface_coeff 0.1
#define initial_radius 1.0
#define c_matrix 0.01
#define c_precip 0.15
#define adjust_avg_c false
#define c_avg 0.004

//initial condition for concentration
template <int dim>
class InitialConditionC : public Function<dim>
{
public:
  double shift;
  InitialConditionC (double _shift) : Function<dim>(1), shift(_shift) {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const Point<dim> &p, const unsigned int component = 0) const
  {

	//return the value of the initial concentration field at point p
	double dx=spanX/( (double)subdivisionsX )/std::pow(2.0,refineFactor);
	double dy=spanY/( (double)subdivisionsY )/std::pow(2.0,refineFactor);
	double dz=spanZ/( (double)subdivisionsZ )/std::pow(2.0,refineFactor);
	double r = 0.0;
    //return 0.02 + 1.0e-3*(2*(0.5 - (double)(std::rand() % 100 )/100.0));
	#if problemDIM==1
	  r=p.operator()(0);
	  return 0.5*(0.12-0.00)*(1-std::tanh((r-spanX/16.0)/(0.1*dx)));
	#elif problemDIM==2
	  // My initial conditions
	  //r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom);
	  //return 0.5*(0.12-0.0)*(1.0-std::tanh((r-spanY/8.0)/(2.0*dy))) +0.02;

	  // Larry's initial conditions (150423f)
	  //r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom);
	  //return 0.5*(0.12-0.0)*(1.0-std::tanh((r-3.0)/(0.5*dy))) +0.01;

	  // Larry's initial conditions (150422p)
	  //r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom);
	  //return 0.5*(0.12-0.0)*(1.0-std::tanh((r-spanX/8.0)/(2.0*dy))) +0.02;

	  // Constant concentration
	  //return avg_Nd;

	  // Ellipsoid
	  r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom
	  		+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom);
	  return 0.5*(c_precip-c_matrix)*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff))) + c_matrix + shift;

	  //double test;
	  //computeIntegral(test);

	  // planar interface
	  //r=(p.operator()(0)-spanX/4.0);
	  //return 0.5*(0.12-0.0)*(1.0-std::tanh((r)/(initial_interface_coeff))) + avg_Nd;

	#elif problemDIM==3
	  // Sphere
	  r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom
		+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom
		+(p.operator()(2)-spanZ/2.0)*(p.operator()(2)-spanZ/2.0)/z_denom);
	  return 0.5*(c_precip-c_matrix)*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff))) + c_matrix + shift;

	  // Constant concentration
	  //return c_matrix + shift;

	  // planar interface
	  //r=sqrt((p.operator()(2))*(p.operator()(2)));
	  //return 0.5*(0.12-avg_Nd)*(1.0-std::tanh((r)/(initial_interface_coeff))) + c_matrix + shift;
	  //return  0.5*(0.12-avg_Nd)*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff))) + c_matrix + shift;

	#endif
  }
};

//initial condition for the structural order parameters
template <int dim>
class InitialConditionN : public Function<dim>
{
public:
  unsigned int index;
  InitialConditionN (const unsigned int _index) : Function<dim>(1), index(_index) {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const Point<dim> &p, const unsigned int component = 0) const
  {
    //set result equal to the structural order paramter initial condition
	  double dx=spanX/( (double)subdivisionsX )/std::pow(2.0,refineFactor);
	  double dy=spanY/( (double)subdivisionsY )/std::pow(2.0,refineFactor);
	  double dz=spanZ/( (double)subdivisionsZ )/std::pow(2.0,refineFactor);
    double r=0.0;
	#if problemDIM==1
    r=p.operator()(0);
	  return 0.5*(1.0-std::tanh((r-spanX/16.0)/(0.1*dx)));
	#elif problemDIM==2
	  if (index==1){
		  // My initial conditions
		  //r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom);
		  //return 0.5*(1.0-std::tanh((r-spanY/8.0)/(2.0*dy)));

		  // Larry's initial conditions (150423f)
		  //r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom);
		  //return 0.5*(1.0-std::tanh((r-3.0)/(0.5*dy)));

		  // Larry's initial conditions (150422p)
//		  r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom);
//		  return 0.5*(1.0-std::tanh((r-spanX/8.0)/(2.0*dy)));

		  r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom
		  		  		+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom);
		  return 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));

		  // planar interface
		  //r=(p.operator()(0)-spanX/4.0);
		  //return 0.5*(1.0-std::tanh((r)/(2.0*dy)));

		  // pure matrix phase
		  //return 0.0;
	  }
	  else if (index==2){
		return 0.0;
	//     double r1=p.distance(Point<dim>(3*spanX/4.0,spanY/4.0));
	//     double r2=p.distance(Point<dim>(spanX/2.0,spanY/2.0));
	//     r=std::min(r1,r2);
	  }
	  else if (index==3){
		return 0.0;
	//     r=p.distance(Point<dim>(spanX/4.0,3*spanY/4.0));
	  }
	  //return 0.5*(1.0-std::tanh((r-spanX/16.0)/(3*dx)));
	#elif problemDIM==3
	  if (index==1){
	  //r=p.distance(Point<dim>(spanX/2.0,spanY/2.0,spanZ/2.0));
	  //return 0.5*(1.0-std::tanh((r-spanX/8.0)/(3*dx)));

		// Sphere
		r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom
		  		+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom
		  		+(p.operator()(2)-spanZ/2.0)*(p.operator()(2)-spanZ/2.0)/z_denom);
		return 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));

		// planar interface
		//r=sqrt((p.operator()(2))*(p.operator()(2)));
		//return 0.5*(1.0-std::tanh((r)/(initial_interface_coeff)));
		//return 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));
	  }
	  else if (index==2){
		return 0.0;
	  }
	  else if (index==3){
		return 0.0;
	  }
	#endif
	return 0.0;
  }
};


//apply initial conditions
template <int dim>
void CoupledCHACMechanicsProblem<dim>::applyInitialConditions()
{


	unsigned int fieldIndex;
	  //call initial condition function for c
	  fieldIndex=this->getFieldIndex("c");
	  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionC<dim>(0.0), *this->solutionSet[fieldIndex]);
	  //call initial condition function for structural order parameters
	  fieldIndex=this->getFieldIndex("n1");
	  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionN<dim>(1), *this->solutionSet[fieldIndex]);
	  if (num_sop > 1){
		  fieldIndex=this->getFieldIndex("n2");
		  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionN<dim>(2), *this->solutionSet[fieldIndex]);
		  if (num_sop > 2){
			  fieldIndex=this->getFieldIndex("n3");
			  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionN<dim>(3), *this->solutionSet[fieldIndex]);
		  }
	  }

	  //set zero intial condition for u
	  fieldIndex=this->getFieldIndex("u");
	  *this->solutionSet[fieldIndex]=0.0;

}

//apply Dirchlet BC function
template <int dim>
void CoupledCHACMechanicsProblem<dim>::applyDirichletBCs(){
  //Set u=0 at all boundaries
  VectorTools::interpolate_boundary_values (*this->dofHandlersSet[this->getFieldIndex("u")],\
					    0, ZeroFunction<dim>(dim), *(ConstraintMatrix*) \
					    this->constraintsSet[this->getFieldIndex("u")]);
}

// Shift the initial concentration so that the average concentration is the desired value
template <int dim>
void CoupledCHACMechanicsProblem<dim>::shiftConcentration()
{
	unsigned int fieldIndex;
	fieldIndex=this->getFieldIndex("c");

	double integrated_concentration;
	computeIntegral(integrated_concentration);

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

	*this->solutionSet[fieldIndex]=0.0;
	VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionC<dim>(shift), *this->solutionSet[fieldIndex]);
	MatrixFreePDE<dim>::solutionSet[fieldIndex]->update_ghost_values();

	computeIntegral(integrated_concentration);
}

//methods to mark boundaries
template <int dim>
void CoupledCHACMechanicsProblem<dim>::markBoundaries(){

	// By default leave all boundaries marked as zero
}
