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

	// Pure virtual method in MatrixFreePDE
	void energyDensity(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list, const dealii::VectorizedArray<double> & JxW_value,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc);

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

	double McV = userInputs.get_model_constant_double(0);
	double Mn1V = userInputs.get_model_constant_double(1);
	dealii::Tensor<2,dim> Kn1 = userInputs.get_model_constant_rank_2_tensor(2);
	double W = userInputs.get_model_constant_double(3);
	bool n_dependent_stiffness = userInputs.get_model_constant_bool(4);
	dealii::Tensor<2,dim> sfts_linear1 = userInputs.get_model_constant_rank_2_tensor(5);
	dealii::Tensor<2,dim> sfts_const1 = userInputs.get_model_constant_rank_2_tensor(6);
	double A2 = userInputs.get_model_constant_double(7);
	double A1 = userInputs.get_model_constant_double(8);
	double A0 = userInputs.get_model_constant_double(9);
	double B2 = userInputs.get_model_constant_double(10);
	double B1 = userInputs.get_model_constant_double(11);
	double B0 = userInputs.get_model_constant_double(12);

	const static unsigned int CIJ_tensor_size =2*dim-1+dim/3;

	// ================================================================

};
