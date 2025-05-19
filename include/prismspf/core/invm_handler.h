// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/types.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

struct VariableAttributes;

/**
 * \brief This class handles the computation and access of the inverted mass matrix for
 * explicit solves.
 */
template <unsigned int dim, unsigned int degree, typename number = double>
class InvmHandler
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;
  using SizeType   = dealii::VectorizedArray<number>;

  /**
   * \brief Constructor.
   */
  explicit InvmHandler(
    const std::map<unsigned int, VariableAttributes> &_variable_attributes);

  /**
   * \brief Initialize.
   */
  void
  initialize(std::shared_ptr<dealii::MatrixFree<dim, number, SizeType>> _data);

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
  const std::map<unsigned int, VariableAttributes> *variable_attributes;

  /**
   * \brief Matrix-free object.
   */
  std::shared_ptr<dealii::MatrixFree<dim, number, SizeType>> data;

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
