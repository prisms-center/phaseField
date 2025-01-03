/*
 * typeDefs.h
 *
 *  Created on: Feb 24, 2017
 *      Author: stephendewitt
 */

// #ifndef INCLUDE_TYPEDEFS_H_
// #define INCLUDE_TYPEDEFS_H_

// define FE system types
#ifndef typeScalar
using typeScalar = dealii::FEEvaluation<dim, degree, degree + 1, 1, double>;
#endif
#ifndef typeVector
using typeVector = dealii::FEEvaluation<dim, degree, degree + 1, dim, double>;
#endif
// define data value types
#ifndef scalarvalueType
using scalarvalueType = dealii::VectorizedArray<double>;
#endif
#ifndef vectorvalueType
using vectorvalueType = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
#endif
#if problemDIM == 1
#  ifndef scalargradType
using scalargradType = dealii::VectorizedArray<double>;
#  endif
#  ifndef vectorgradType
using vectorgradType = dealii::VectorizedArray<double>;
#  endif
#  ifndef vectorhessType
using vectorhessType = dealii::VectorizedArray<double>;
#  endif
#else
#  ifndef scalargradType
using scalargradType = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
#  endif
#  ifndef scalarhessType
using scalarhessType = dealii::Tensor<2, dim, dealii::VectorizedArray<double>>;
#  endif
#  ifndef vectorgradType
using vectorgradType = dealii::Tensor<2, dim, dealii::VectorizedArray<double>>;
#  endif
#  ifndef vectorhessType
using vectorhessType = dealii::Tensor<3, dim, dealii::VectorizedArray<double>>;
#  endif
#endif

// #endif /* INCLUDE_TYPEDEFS_H_ */
