#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:

	// Constructor, which calls the MatrixFreePDE constructor
	customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs), userInputs(_userInputs) {};

private:
	#include "../../include/typeDefs.h"

	const userInputParameters<dim> userInputs;

	// Pure virtual method in MatrixFreePDE
	void residualRHS(const std::vector<modelVariable<dim> > & modelVarList,
			  	  	 std::vector<modelResidual<dim> > & modelResidualsList,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Pure virtual method in MatrixFreePDE
	void residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
	  		  	  	 modelResidual<dim> & modelRes,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Pure virtual method in MatrixFreePDE
	void energyDensity(const std::vector<modelVariable<dim> > & modelVarList, const dealii::VectorizedArray<double> & JxW_value,
			  	  	 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc);

	// Virtual method in MatrixFreePDE that we override if we need postprocessing (and a postprocess.h file exists)
	#ifdef POSTPROCESS_FILE_EXISTS
	void postProcessedFields(const std::vector<modelVariable<dim> > & modelVariablesList,
				 	std::vector<modelResidual<dim> > & modelResidualsList,
				 	const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
	#endif

	// Virtual method in MatrixFreePDE that we override if we need nucleation (and a nucleation.h file exists)
	#ifdef NUCLEATION_FILE_EXISTS
	double getNucleationProbability(variableValueContainer variable_value, double dV) const;
	#endif

	// ================================================================
	// Methods specific to this subclass
	// ================================================================


	// ================================================================
	// Model constants specific to this subclass
	// ================================================================

	double McV = userInputs.get_model_constant_double(0);
	double MnV = userInputs.get_model_constant_double(1);
	double KnV = userInputs.get_model_constant_double(2);
	double c_avg = userInputs.get_model_constant_double(3);
	double W_barrier = userInputs.get_model_constant_double(4);
	double A0 = userInputs.get_model_constant_double(5);
	double A2 = userInputs.get_model_constant_double(6);
	double calmin = userInputs.get_model_constant_double(7);
	double B0 = userInputs.get_model_constant_double(8);
	double B2 = userInputs.get_model_constant_double(9);
	double cbtmin = userInputs.get_model_constant_double(10);
	double k1 = userInputs.get_model_constant_double(11);
	double k2 = userInputs.get_model_constant_double(12);
	double tau = userInputs.get_model_constant_double(13);
	double epsilon = userInputs.get_model_constant_double(14);

	double interface_coeff=std::sqrt(2.0*KnV/W_barrier);

	// ================================================================

};
