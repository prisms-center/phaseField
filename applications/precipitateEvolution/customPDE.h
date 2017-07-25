#include "../../include/matrixFreePDE.h"
//#include "../../include/variableAttributeLoader.h"

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

	double McV = userInputs.get_model_constant_double("McV");
	double Mn1V = userInputs.get_model_constant_double("Mn1V");
	double Mn2V = userInputs.get_model_constant_double("Mn2V");
	double Mn3V = userInputs.get_model_constant_double("Mn3V");
	dealii::Tensor<2,dim> Kn1 = userInputs.get_model_constant_rank_2_tensor("Kn1");
	dealii::Tensor<2,dim> Kn2 = userInputs.get_model_constant_rank_2_tensor("Kn2");
	dealii::Tensor<2,dim> Kn3 = userInputs.get_model_constant_rank_2_tensor("Kn3");
	bool n_dependent_stiffness = userInputs.get_model_constant_bool("n_dependent_stiffness");
	dealii::Tensor<2,dim> sfts_linear1 = userInputs.get_model_constant_rank_2_tensor("sfts_linear1");
	dealii::Tensor<2,dim> sfts_const1 = userInputs.get_model_constant_rank_2_tensor("sfts_const1");
	dealii::Tensor<2,dim> sfts_linear2 = userInputs.get_model_constant_rank_2_tensor("sfts_linear2");
	dealii::Tensor<2,dim> sfts_const2 = userInputs.get_model_constant_rank_2_tensor("sfts_const2");
	dealii::Tensor<2,dim> sfts_linear3 = userInputs.get_model_constant_rank_2_tensor("sfts_linear3");
	dealii::Tensor<2,dim> sfts_const3 = userInputs.get_model_constant_rank_2_tensor("sfts_const3");
	double A4 = userInputs.get_model_constant_double("A4");
	double A3 = userInputs.get_model_constant_double("A3");
	double A2 = userInputs.get_model_constant_double("A2");
	double A1 = userInputs.get_model_constant_double("A1");
	double A0 = userInputs.get_model_constant_double("A0");
	double B2 = userInputs.get_model_constant_double("B2");
	double B1 = userInputs.get_model_constant_double("B1");
	double B0 = userInputs.get_model_constant_double("B0");

	const static unsigned int CIJ_tensor_size =2*dim-1+dim/3;
	dealii::Tensor<2,CIJ_tensor_size> CIJ_Mg = userInputs.get_model_constant_elasticity_tensor("CIJ_Mg");
	dealii::Tensor<2,CIJ_tensor_size> CIJ_Beta = userInputs.get_model_constant_elasticity_tensor("CIJ_Beta");



	// ================================================================

};
