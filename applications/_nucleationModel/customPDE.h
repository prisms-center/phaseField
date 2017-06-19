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

	// Virtual method in MatrixFreePDE we choose to override
	void getNucleiList ();

	const userInputParameters<dim> userInputs;

	// ================================================================
	// Model constants
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

	// ----------------------------------------------------------------
	// Nucleation methods specific to this subclass
	// ----------------------------------------------------------------

	// Contains nucleation probability that varies between applications, no MatrixFreePDE member access
	double getNucleationProbability(variableValueContainer variable_value, double dV) const;


};

#endif /* APPLICATIONS_ALLENCAHN_CUSTOMPDE_H_ */
