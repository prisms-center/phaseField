/*
 * typeDefs.h
 *
 *  Created on: Feb 24, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_TYPEDEFS_H_
#define INCLUDE_TYPEDEFS_H_

//define data types
#ifndef scalarType
typedef dealii::VectorizedArray<double> scalarType;
#endif
#ifndef vectorType
typedef dealii::parallel::distributed::Vector<double> vectorType;
#endif
//define FE system types
#ifndef typeScalar
typedef dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,1,double>           typeScalar;
#endif
#ifndef typeVector
typedef dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,problemDIM,double>  typeVector;
#endif
//define data value types
#ifndef scalarvalueType
typedef dealii::VectorizedArray<double> scalarvalueType;
#endif
#ifndef vectorvalueType
typedef dealii::Tensor<1, problemDIM, dealii::VectorizedArray<double> > vectorvalueType;
#endif
#if problemDIM==1
#ifndef scalargradType
typedef dealii::VectorizedArray<double> scalargradType;
#endif
#ifndef vectorgradType
typedef dealii::VectorizedArray<double> vectorgradType;
#endif
#ifndef vectorhessType
typedef dealii::VectorizedArray<double> vectorhessType;
#endif
#else
#ifndef scalargradType
typedef dealii::Tensor<1, problemDIM, dealii::VectorizedArray<double> > scalargradType;
#endif
#ifndef scalarhessType
typedef dealii::Tensor<2,problemDIM,dealii::VectorizedArray<double> > scalarhessType;
#endif
#ifndef vectorgradType
typedef dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > vectorgradType;
#endif
#ifndef vectorhessType
typedef dealii::Tensor<3, problemDIM, dealii::VectorizedArray<double> > vectorhessType;
#endif
#endif

#endif /* INCLUDE_TYPEDEFS_H_ */
