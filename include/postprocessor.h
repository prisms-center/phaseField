/*
 * postprocessor.h
 *
 *  Created on: Mar 29, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_POSTPROCESSOR_H_
#define INCLUDE_POSTPROCESSOR_H_

#include "dealIIheaders.h"
//#include "typeDefs.h"
#include "matrixFreePDE.h"
#include "postProcessedFields.h"

template <int dim, int degree>
class PostProcessor
{
	public:
	PostProcessor(userInputParameters);
	void computePostProcessedFields(dealii::MatrixFree<dim,double> matrixFreeObject,
							const std::vector<vectorType*> &solutionSet,
							std::vector<vectorType*> &postProcessedSet);
	void getPostProcessedFields(const dealii::MatrixFree<dim,double> &data,
							std::vector<vectorType*> &dst,
							const std::vector<vectorType*> &src,
							const std::pair<unsigned int,unsigned int> &cell_range);
//	void postProcessedFields(const std::vector<modelVariable<dim> > & modelVariablesList,
//							std::vector<modelResidual<dim> > & modelResidualsList,
//							const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
	private:
	userInputParameters userInputs;

};



#endif /* INCLUDE_POSTPROCESSOR_H_ */
