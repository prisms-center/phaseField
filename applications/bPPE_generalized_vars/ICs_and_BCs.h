//initial condition for concentration
template <int dim>
class InitialConditionC : public Function<dim>
{
public:
  InitialConditionC () : Function<dim>(1) {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const Point<dim> &p, const unsigned int component = 0) const
  {
    //return the value of the initial concentration field at point p
    return 0.04; //+ 1.0e-3*(2*(0.5 - (double)(std::rand() % 100 )/100.0));
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
    double dx=spanX/((double) subdivisionsX)/std::pow(2.0,refineFactor);
    double r=0.0;
#if problemDIM==2
    if (index==1){
      double r1=p.distance(Point<dim>(spanX/3.0,spanY/3.0));
      double r2=p.distance(Point<dim>(2*spanX/3.0,2*spanY/3.0));
      r=std::min(r1,r2);
    }
    else if (index==2){
      r=p.distance(Point<dim>(3*spanX/4.0,spanY/4.0));
    }
    else if (index==3){
      r=p.distance(Point<dim>(spanX/4.0,3*spanY/4.0));
    }
    return 0.5*(1.0-std::tanh((r-spanX/16.0)/(dx)));
#endif
  }
};


//apply initial conditions
template <int dim>
void CoupledCHACMechanicsProblem<dim>::applyInitialConditions()
{
  unsigned int fieldIndex;
  //call initial condition function for c
  fieldIndex=this->getFieldIndex("c");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionC<dim>(), *this->solutionSet[fieldIndex]);
  //call initial condition function for structural order parameters
  fieldIndex=this->getFieldIndex("n1");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionN<dim>(1), *this->solutionSet[fieldIndex]);
  fieldIndex=this->getFieldIndex("n2");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionN<dim>(2), *this->solutionSet[fieldIndex]);
  fieldIndex=this->getFieldIndex("n3");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionN<dim>(3), *this->solutionSet[fieldIndex]);
  //set zero initial condition for u
  fieldIndex=this->getFieldIndex("u");
  *this->solutionSet[fieldIndex]=0.0;
}

template <int dim>
class vectorBCFunction : public Function<dim>
{
public:
	vectorBCFunction(std::vector<double> BC_values);
  void vector_value (const Point<dim> &p, Vector<double>   &values) const
  {
	  for (unsigned int component=0; component < dim; component++){
		  values(component) = BC_values[component];
	  }
  }
private:
  std::vector<double> BC_values;
};

template <int dim>
vectorBCFunction<dim>::vectorBCFunction(std::vector<double> input_values){
	BC_values = input_values;
}


//apply Dirchlet BC function
template <int dim>
void CoupledCHACMechanicsProblem<dim>::applyDirichletBCs(){
  // First, get the variable index of the current field
  unsigned int var_index;
  unsigned int field_number = 0;
  for (unsigned int i=0; i<num_var; i++){

	  if (field_number == this->currentFieldIndex){
		  var_index = i;
	  }

	  if (var_type[var_index] == "SCALAR"){
		  field_number++;
	  }
	  else {
		  field_number+=dim;
	  }
  }

  if (var_type[var_index] == "SCALAR"){
	  for (unsigned int direction = 0; direction < 2*dim; direction++){
		  if (BC_list[this->currentFieldIndex].var_BC_type[direction] == "DIRICHLET"){
			  VectorTools::interpolate_boundary_values (*this->dofHandlersSet[this->currentFieldIndex],\
					  direction, ConstantFunction<dim>(BC_list[this->currentFieldIndex].var_BC_val[direction]), *(ConstraintMatrix*) \
					  this->constraintsSet[this->currentFieldIndex]);
		  }
	  }
  }
  else {
	  for (unsigned int direction = 0; direction < 2*dim; direction++){

		  std::vector<double> BC_values;
		  for (unsigned int component=0; component < dim; component++){
			  BC_values.push_back(BC_list[this->currentFieldIndex+component].var_BC_val[direction]);
		  }


		  if (BC_list[this->currentFieldIndex].var_BC_type[direction] == "DIRICHLET"){
			  VectorTools::interpolate_boundary_values (*this->dofHandlersSet[this->currentFieldIndex],\
							direction, vectorBCFunction<dim>(BC_values), *(ConstraintMatrix*) \
							this->constraintsSet[this->currentFieldIndex]);
		  }

	  }
  }
}

//methods to mark boundaries
template <int dim>
void CoupledCHACMechanicsProblem<dim>::markBoundaries(int field_index){

	std::vector<double> domain_size;
	domain_size.push_back(spanX);
	domain_size.push_back(spanY);
	domain_size.push_back(spanZ);

	typename Triangulation<dim>::cell_iterator
	cell = MatrixFreePDE<dim>::triangulation.begin (),
	endc = MatrixFreePDE<dim>::triangulation.end();

	for (; cell!=endc; ++cell){

		// Mark all of the faces
		for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell;++face_number){
			for (unsigned int i=0; i<dim; i++){
				if ( std::fabs(cell->face(face_number)->center()(i) - (0)) < 1e-12 ){
					cell->face(face_number)->set_boundary_indicator (2*i);
				}
				else if (std::fabs(cell->face(face_number)->center()(i) - (domain_size[i])) < 1e-12){
					cell->face(face_number)->set_boundary_indicator (2*i+1);
				}

			}
		}
	}
}

// Input the boundary conditions
template <int dim>
void CoupledCHACMechanicsProblem<dim>::setBCs(){
	// Set the BCs for the problem variables in numerical order
	// Four input arguments set the same BC on the entire boundary
	// Two plus two times the number of dimensions inputs sets separate BCs on each face of the domain

	inputBCs(0,0,"ZERO_DERIVATIVE",0);

	inputBCs(1,0,"ZERO_DERIVATIVE",0);

	inputBCs(2,0,"ZERO_DERIVATIVE",0);

	inputBCs(3,0,"ZERO_DERIVATIVE",0);

	inputBCs(4,0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0);
	inputBCs(4,1,"DIRICHLET",0.0);

}

// Input the boundary conditions
template <int dim>
void CoupledCHACMechanicsProblem<dim>::inputBCs(int var, int component, std::string BC_type_dim1_min, double BC_value_dim1_min,
		std::string BC_type_dim1_max, double BC_value_dim1_max, std::string BC_type_dim2_min, double BC_value_dim2_min,
		std::string BC_type_dim2_max, double BC_value_dim2_max){

	varBCs<dim> newBC;
	newBC.var_BC_type.push_back(BC_type_dim1_min);
	newBC.var_BC_type.push_back(BC_type_dim1_max);
	newBC.var_BC_type.push_back(BC_type_dim2_min);
	newBC.var_BC_type.push_back(BC_type_dim2_max);

	newBC.var_BC_val.push_back(BC_value_dim1_min);
	newBC.var_BC_val.push_back(BC_value_dim1_max);
	newBC.var_BC_val.push_back(BC_value_dim2_min);
	newBC.var_BC_val.push_back(BC_value_dim2_max);

	BC_list.push_back(newBC);
}

template <int dim>
void CoupledCHACMechanicsProblem<dim>::inputBCs(int var, int component, std::string BC_type, double BC_value){

	varBCs<dim> newBC;
	for (unsigned int face=0; face<dim*2; face++){
		newBC.var_BC_type.push_back(BC_type);
		newBC.var_BC_val.push_back(BC_value);
	}

	BC_list.push_back(newBC);
}
