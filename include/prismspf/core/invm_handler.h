// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef invm_handler_h
#define invm_handler_h

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/config.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <map>
#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class handles the computation and access of the inverted mass matrix for
 * explicit solves.
 */
template <int dim, int degree>
class invmHandler
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<double>;
  using value_type = double;
  using size_type  = dealii::VectorizedArray<double>;

  /**
   * \brief Constructor.
   */
  explicit invmHandler(
    const std::map<unsigned int, variableAttributes> &_variable_attributes);

  /**
   * \brief Destructor.
   */
  ~invmHandler() = default;

  /**
   * \brief Initialize.
   */
  void
  initialize(std::shared_ptr<dealii::MatrixFree<dim, double, size_type>> _data);

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
   * \brief  Getter function for the mass matrix for the given field index (constant
   * reference).
   */
  [[nodiscard]] const VectorType &
  get_invm(const unsigned int &index) const;

private:
  /**
   * \brief Variable attributes. This is used to determine the proper return type for the
   * invm when given a field index.
   */
  const std::map<unsigned int, variableAttributes> &variable_attributes;

  /**
   * \brief Matrix-free object.
   */
  std::shared_ptr<dealii::MatrixFree<dim, double, size_type>> data;

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
  unsigned int scalar_index;

  /**
   * \brief Field index of the first occuring vector field. This is the index for which we
   * attached the FEEvaluation objects to evaluate and initialize the invm vector.
   */
  unsigned int vector_index;
};

PRISMS_PF_END_NAMESPACE

#endif