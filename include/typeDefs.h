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
using scalarType = dealii::VectorizedArray<double>;
#endif
#ifndef vectorType
using vectorType = dealii::LinearAlgebra::distributed::Vector<double>;
#endif
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
