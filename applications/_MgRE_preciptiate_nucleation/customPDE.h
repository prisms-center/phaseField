/*
 * customPDE.h
 *
 *  Created on: Feb 24, 2017
 *      Author: stephendewitt
 */

#ifndef APPLICATIONS_ALLENCAHN_CUSTOMPDE_H_
#define APPLICATIONS_ALLENCAHN_CUSTOMPDE_H_

#include "../../include/matrixFreePDE.h"
#include "../../include/parallelNucleationList.h"
#include "../../include/nucleus.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:

	// Constructor, which calls the MatrixFreePDE constructor
	customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs), userInputs(_userInputs) {};

private:
	#include "../../include/typeDefs.h"

	const static unsigned int CIJ_tensor_size =2*dim-1+dim/3;

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

	// Virtual method in MatrixFreePDE that we override
	void postProcessedFields(const std::vector<modelVariable<dim> > & modelVariablesList,
				 	std::vector<modelResidual<dim> > & modelResidualsList,
				 	const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	const userInputParameters<dim> userInputs;

	// ================================================================
	// Model constants
	// ================================================================

	double McV = userInputs.get_model_constant_double(0);
	double Mn1V = userInputs.get_model_constant_double(1);
	double Mn2V = userInputs.get_model_constant_double(2);
	double Mn3V = userInputs.get_model_constant_double(3);
	dealii::Tensor<2,dim> Kn1 = userInputs.get_model_constant_rank_2_tensor(4);
	dealii::Tensor<2,dim> Kn2 = userInputs.get_model_constant_rank_2_tensor(5);
	dealii::Tensor<2,dim> Kn3 = userInputs.get_model_constant_rank_2_tensor(6);
	double W = userInputs.get_model_constant_double(7);
	bool n_dependent_stiffness = userInputs.get_model_constant_bool(8);
	dealii::Tensor<2,dim> sfts_linear1 = userInputs.get_model_constant_rank_2_tensor(9);
	dealii::Tensor<2,dim> sfts_const1 = userInputs.get_model_constant_rank_2_tensor(10);
	dealii::Tensor<2,dim> sfts_linear2 = userInputs.get_model_constant_rank_2_tensor(11);
	dealii::Tensor<2,dim> sfts_const2 = userInputs.get_model_constant_rank_2_tensor(12);
	dealii::Tensor<2,dim> sfts_linear3 = userInputs.get_model_constant_rank_2_tensor(13);
	dealii::Tensor<2,dim> sfts_const3 = userInputs.get_model_constant_rank_2_tensor(14);
	double A2 = userInputs.get_model_constant_double(15);
	double A1 = userInputs.get_model_constant_double(16);
	double A0 = userInputs.get_model_constant_double(17);
	double B2 = userInputs.get_model_constant_double(18);
	double B1 = userInputs.get_model_constant_double(19);
	double B0 = userInputs.get_model_constant_double(20);

	// ================================================================

	// ----------------------------------------------------------------
	// Nucleation methods specific to this subclass
	// ----------------------------------------------------------------

	// Contains nucleation probability that varies between applications, no MatrixFreePDE member access
	double getNucleationProbability(variableValueContainer variable_value, double dV) const;

	void seedNucleus(const dealii::Point<dim, dealii::VectorizedArray<double> > & q_point_loc,
						std::vector<dealii::VectorizedArray<double> > & source_terms,
						dealii::VectorizedArray<double> & gamma) const;

};

#endif /* APPLICATIONS_ALLENCAHN_CUSTOMPDE_H_ */
