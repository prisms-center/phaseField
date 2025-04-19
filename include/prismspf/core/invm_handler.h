// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/types.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
class userInputParameters;

struct variableAttributes;

/**
 * \brief This class handles the computation and access of the inverted mass matrix for
 * explicit solves.
 */
template <int dim, int degree, typename number = double>
class invmHandler
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;
  using size_type  = dealii::VectorizedArray<number>;

  /**
   * \brief Constructor.
   */
  explicit invmHandler(
    const std::map<unsigned int, variableAttributes> &_variable_attributes);

  /**
   * \brief Initialize.
   */
  void
  initialize(std::shared_ptr<dealii::MatrixFree<dim, number, size_type>> _data);

  /**
   * \brief Compute the mass matrix for scalar/vector fields.
   */
  void
  compute_invm();

  /**
   * \brief Recompute the mass matrix for scalar/vector fields. This just points to
   * compute_invm() and is used for style.
   */
  void
  recompute_invm();

  /**
   * \brief Getter function for the mass matrix for the given field index (constant
   * reference).
   */
  [[nodiscard]] const VectorType &
  get_invm(const unsigned int &index) const;

  /**
   * \brief Clear the data so we have something that resembles the base constructor.
   */
  void
  clear();

private:
  /**
   * \brief Compute the invm for scalar fields.
   */
  void
  compute_scalar_invm();

  /**
   * \brief Compute the invm for vector fields.
   */
  void
  compute_vector_invm();

  /**
   * \brief Variable attributes. This is used to determine the proper return type for the
   * invm when given a field index.
   */
  const std::map<unsigned int, variableAttributes> *variable_attributes;

  /**
   * \brief Matrix-free object.
   */
  std::shared_ptr<dealii::MatrixFree<dim, number, size_type>> data;

  /**
   * \brief Inverse of the mass matrix for scalar fields.
   */
  VectorType invm_scalar;

  /**
   * \brief Inverse of the mass matrix for vector fields.
   */
  VectorType invm_vector;

  /**
   * \brief Whether a scalar invm is needed.
   */
  bool scalar_needed = false;

  /**
   * \brief Whether a vector invm is needed.
   */
  bool vector_needed = false;

  /**
   * \brief Field index of the first occuring scalar field. This is the index for which we
   * attached the FEEvaluation objects to evaluate and initialize the invm vector.
   */
  unsigned int scalar_index = numbers::invalid_index;

  /**
   * \brief Field index of the first occuring vector field. This is the index for which we
   * attached the FEEvaluation objects to evaluate and initialize the invm vector.
   */
  unsigned int vector_index = numbers::invalid_index;

  /**
   * \brief Tolerance for minimum value of the mass matrix when inverting.
   */
  number tolerance = defaults::mesh_tolerance;
};

PRISMS_PF_END_NAMESPACE
