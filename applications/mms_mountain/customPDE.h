#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
	customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {};

private:
	#include "../../include/typeDefs.h"

	const userInputParameters<dim> userInputs;

	// Pure virtual method in MatrixFreePDE
	void residualRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Pure virtual method in MatrixFreePDE
	void residualLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Virtual method in MatrixFreePDE that we override if we need postprocessing
	#ifdef POSTPROCESS_FILE_EXISTS
	void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
					const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
	#endif

	// Virtual method in MatrixFreePDE that we override if we need nucleation
	#ifdef NUCLEATION_FILE_EXISTS
	double getNucleationProbability(variableValueContainer variable_value, double dV) const;
	#endif

	// ================================================================
	// Methods specific to this subclass
	// ================================================================


	// ================================================================
	// Model constants specific to this subclass
	// ================================================================

	double MnV = userInputs.get_model_constant_double("MnV");
	double KnV = userInputs.get_model_constant_double("KnV");

	double A1 = userInputs.get_model_constant_double("A1");
	double B1 = userInputs.get_model_constant_double("B1");
	double C1 = userInputs.get_model_constant_double("C1");
	double A2 = userInputs.get_model_constant_double("A2");
	double B2 = userInputs.get_model_constant_double("B2");
	double D2 = userInputs.get_model_constant_double("D2");

	double r0 = userInputs.get_model_constant_double("r0");

	double delta = 1.0/std::sqrt(KnV*2.0);


	// ================================================================

};
