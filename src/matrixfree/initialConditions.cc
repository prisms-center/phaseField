//methods to apply initial conditions 

#include "../../include/matrixFreePDE.h"
#include "../../include/IntegrationTools/PField.hh"

#ifdef enablePFields
#if enablePFields == true
template <int dim>
class InitialConditionPField : public Function<dim>
{
public:
  unsigned int index;
  Vector<double> values;
  typedef PRISMS::PField<double*, double, 2> ScalarField2D;
  ScalarField2D &inputField;

  InitialConditionPField (const unsigned int _index, ScalarField2D &_inputField) : Function<dim>(1), index(_index), inputField(_inputField) {}

  double value (const Point<dim> &p, const unsigned int component = 0) const
  {
	  double scalar_IC;

	  double coord[dim];
	  for (unsigned int i = 0; i < dim; i++){
		  coord[i] = p(i);
	  }

	  scalar_IC = inputField(coord);

	  return scalar_IC;
  }
};

#endif
#else
#define enablePFields false
#endif

//methods to apply initial conditions
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::applyInitialConditions(){

for (unsigned int var_index=0; var_index < userInputs.number_of_variables; var_index++){
	if (userInputs.load_ICs[var_index] == false){
		if (userInputs.var_type[var_index] == "SCALAR"){
			VectorTools::interpolate (*dofHandlersSet[var_index], InitialCondition<dim>(var_index), *solutionSet[var_index]);
		}
		else {
			VectorTools::interpolate (*dofHandlersSet[var_index], InitialConditionVec<dim>(var_index), *solutionSet[var_index]);
		}
	}
	else{
		#if enablePFields == true
		// Declare the PField types and containers
		typedef PRISMS::PField<double*, double, 2> ScalarField2D;
		typedef PRISMS::Body<double*, 2> Body2D;
		Body2D body;

		// Create the filename of the the file to be loaded
		std::string filename;
		if (userInputs.load_serial_file[var_index] == true){
			filename = userInputs.load_file_name[var_index] + ".vtk";
		}
		else {
			int proc_num = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
			std::ostringstream conversion;
			conversion << proc_num;
			filename = userInputs.load_file_name[var_index] + "." + conversion.str() + ".vtk";
		}

		// Load the data from the file using a PField
		body.read_vtk(filename);
		ScalarField2D &conc = body.find_scalar_field(userInputs.load_field_name[var_index]);
		if (userInputs.var_type[var_index] == "SCALAR"){
			VectorTools::interpolate (*dofHandlersSet[var_index], InitialConditionPField<dim>(var_index,conc), *solutionSet[var_index]);
		}
		else {
			std::cout << "PRISMS-PF Error: Cannot load vector fields. Loading initial conditions from file is currently limited to scalar fields" << std::endl;
		}
		#else
		std::cout << "PRISMS-PF Error: The parameter \"enablePFields\" must be set to true to load initial conditions from file." << std::endl;
		#endif
	}
}
}




// =================================================================================

// I don't think vector fields are implemented in PFields yet
//template <int dim>
//class InitialConditionPFieldVec : public Function<dim>
//{
//public:
//  unsigned int index;
//  Vector<double> values;
//  typedef PRISMS::PField<double*, double, 2> ScalarField2D;
//  ScalarField2D &inputField;
//
//  InitialConditionPFieldVec (const unsigned int _index, ScalarField2D &_inputField) : Function<dim>(1), index(_index), inputField(_inputField) {}
//
//  void vector_value (const Point<dim> &p,Vector<double> &vector_IC) const
//  {
//	  double coord[dim];
//	  for (unsigned int i = 0; i < dim; i++){
//		  coord[i] = p(i);
//	  }
//
//	  vector_IC = inputField(coord);
//  }
//};

#include "../../include/matrixFreePDE_template_instantiations.h"
