/*
 * customPDE.h
 *
 *  Created on: Feb 24, 2017
 *      Author: stephendewitt
 */

#ifndef APPLICATIONS_PRECIPITATEEVOLUTION_CUSTOMPDE_H_
#define APPLICATIONS_PRECIPITATEEVOLUTION_CUSTOMPDE_H_

#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
	customPDE(userInputParameters _userInputs): MatrixFreePDE<dim,degree>(_userInputs) {};

	void setBCs();

private:

	const static unsigned int CIJ_tensor_size =2*dim-1+dim/3;

	void residualRHS(const std::vector<modelVariable<dim> > & modelVarList,
			  	  	 std::vector<modelResidual<dim> > & modelResidualsList,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	void residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
	  		  	  	 modelResidual<dim> & modelRes,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	void energyDensity(const std::vector<modelVariable<dim> > & modelVarList, const dealii::VectorizedArray<double> & JxW_value,
			  	  	 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc);

};



#endif /* APPLICATIONS_PRECIPITATEEVOLUTION_CUSTOMPDE_H_ */
