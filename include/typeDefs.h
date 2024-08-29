/*
 * typeDefs.h
 *
 *  Created on: Feb 24, 2017
 *      Author: stephendewitt
 */

// #ifndef INCLUDE_TYPEDEFS_H_
// #define INCLUDE_TYPEDEFS_H_

// #include <deal.II/base/quadrature.h>
// #include <deal.II/base/timer.h>
// #include <deal.II/lac/vector.h>
// #include <deal.II/lac/affine_constraints.h>
// #include <deal.II/fe/fe_system.h>
// #include <deal.II/fe/fe_q.h>
// #include <deal.II/fe/fe_values.h>
// #include <deal.II/grid/tria.h>
// #include <deal.II/grid/tria_accessor.h>
// #include <deal.II/grid/tria_iterator.h>
// #include <deal.II/grid/grid_tools.h>
// #include <deal.II/dofs/dof_tools.h>
// #include <deal.II/dofs/dof_handler.h>
//  #include <deal.II/numerics/vector_tools.h>
//  #include <deal.II/lac/la_parallel_vector.h>
//  #include <deal.II/matrix_free/matrix_free.h>
//  #include <deal.II/matrix_free/fe_evaluation.h>
//  #include <deal.II/base/config.h>
//  #include <deal.II/base/exceptions.h>
//  #include <deal.II/distributed/tria.h>
// #include <deal.II/distributed/solution_transfer.h>
// #include <deal.II/grid/manifold_lib.h>

// define data types
#ifndef scalarType
typedef dealii::VectorizedArray<double> scalarType;
#endif
#ifndef vectorType
typedef dealii::LinearAlgebra::distributed::Vector<double> vectorType;
#endif
// define FE system types
#ifndef typeScalar
typedef dealii::FEEvaluation<dim, degree, degree + 1, 1, double> typeScalar;
#endif
#ifndef typeVector
typedef dealii::FEEvaluation<dim, degree, degree + 1, dim, double> typeVector;
#endif
// define data value types
#ifndef scalarvalueType
typedef dealii::VectorizedArray<double> scalarvalueType;
#endif
#ifndef vectorvalueType
typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double>> vectorvalueType;
#endif
#if problemDIM == 1
#  ifndef scalargradType
typedef dealii::VectorizedArray<double> scalargradType;
#  endif
#  ifndef vectorgradType
typedef dealii::VectorizedArray<double> vectorgradType;
#  endif
#  ifndef vectorhessType
typedef dealii::VectorizedArray<double> vectorhessType;
#  endif
#else
#  ifndef scalargradType
typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double>> scalargradType;
#  endif
#  ifndef scalarhessType
typedef dealii::Tensor<2, dim, dealii::VectorizedArray<double>> scalarhessType;
#  endif
#  ifndef vectorgradType
typedef dealii::Tensor<2, dim, dealii::VectorizedArray<double>> vectorgradType;
#  endif
#  ifndef vectorhessType
typedef dealii::Tensor<3, dim, dealii::VectorizedArray<double>> vectorhessType;
#  endif
#endif

// #endif /* INCLUDE_TYPEDEFS_H_ */
