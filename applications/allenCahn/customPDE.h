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
	// void residualRHS(const std::vector<modelVariable<dim> > & modelVarList,
	// 				 std::vector<modelResidual<dim> > & modelResidualsList,
	// 				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
	void residualRHS(const variableContainer<dim,dealii::VectorizedArray<double> > & variable_list,
					 std::vector<modelResidual<dim> > & modelResidualsList,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Pure virtual method in MatrixFreePDE
	void residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
					 modelResidual<dim> & modelRes,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Pure virtual method in MatrixFreePDE
	void energyDensity(const std::vector<modelVariable<dim> > & modelVarList, const dealii::VectorizedArray<double> & JxW_value,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc);

	// Virtual method in MatrixFreePDE that we override if we need postprocessing
	#ifdef POSTPROCESS_FILE_EXISTS
	void postProcessedFields(const std::vector<modelVariable<dim> > & modelVariablesList,
					std::vector<modelResidual<dim> > & modelResidualsList,
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

	double MnV = userInputs.get_model_constant_double(0);
	double KnV = userInputs.get_model_constant_double(1);

	// ================================================================

};
