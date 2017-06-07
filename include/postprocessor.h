/*
 * postprocessor.h
 *
 *  Created on: Mar 29, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_POSTPROCESSOR_H_
#define INCLUDE_POSTPROCESSOR_H_

#include "dealIIheaders.h"
#include "varBCs.h"
#include "userInputParameters.h"

////define data types
#ifndef scalarType
typedef dealii::VectorizedArray<double> scalarType;
#endif
#ifndef vectorType
typedef dealii::parallel::distributed::Vector<double> vectorType;
#endif

#include "postProcessedFields.h"

template <int dim, int degree>
class PostProcessor
{
	public:
	PostProcessor(userInputParameters<dim>);
	void computePostProcessedFields(dealii::MatrixFree<dim,double> matrixFreeObject,
							const std::vector<vectorType*> &solutionSet,
							std::vector<vectorType*> &postProcessedSet);
	void getPostProcessedFields(const dealii::MatrixFree<dim,double> &data,
							std::vector<vectorType*> &dst,
							const std::vector<vectorType*> &src,
							const std::pair<unsigned int,unsigned int> &cell_range);
	private:
	userInputParameters<dim> userInputs;

};



#endif /* INCLUDE_POSTPROCESSOR_H_ */
