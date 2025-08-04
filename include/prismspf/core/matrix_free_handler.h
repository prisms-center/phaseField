// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

template <unsigned int dim>
class MGInfo;

template <unsigned int dim, unsigned int degree, typename number>
class ConstraintHandler;

template <unsigned int dim>
class DofHandler;

/**
 * @brief This class is a simple wrapper around deal.II's MatrixFree class reduced in a
 * way for PRISMS-PF.
 */
template <unsigned int dim, typename number>
class MatrixFreeHandler
{
public:
  /**
   * @brief Constructor.
   */
  MatrixFreeHandler();

  /**
   * @brief Reinitialize the matrix-free object with the same quad rule.
   */
  void
  reinit(const dealii::Mapping<dim>              &mapping,
         const dealii::DoFHandler<dim>           &dof_handler,
         const dealii::AffineConstraints<number> &constraint,
         const dealii::Quadrature<1>             &quad);

  /**
   * @brief Reinitialize the matrix-free object with the same quad rule.
   */
  void
  reinit(const dealii::Mapping<dim>                                   &mapping,
         const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
         const std::vector<const dealii::AffineConstraints<number> *> &constraint,
         const dealii::Quadrature<1>                                  &quad);

  /**
   * @brief Reinitialize the matrix-free object with the different quad rule.
   */
  void
  reinit(const dealii::Mapping<dim>                                   &mapping,
         const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
         const std::vector<const dealii::AffineConstraints<number> *> &constraint,
         const std::vector<dealii::Quadrature<1>>                     &quad);

  /**
   * @brief Getter function for the matrix-free object (shared ptr).
   */
  [[nodiscard]] std::shared_ptr<
    dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>
  get_matrix_free() const;

private:
  /**
   * @brief Matrix-free object that collects data to be used in cell loop operations.
   */
  std::shared_ptr<dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>
    matrix_free_object;

  /**
   * @brief Additional data scheme
   */
  typename dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>::
    AdditionalData additional_data;

  /**
   * @brief Whether the matrix-free object has been initialized.
   */
  bool is_initialized = false;
};

/**
 * @brief This class is a simple container for the matrix-free objects used in PRISMS-PF.
 */
template <unsigned int dim, typename number>
class MatrixFreeContainer
{
public:
  /**
   * @brief Constructor.
   */
  explicit MatrixFreeContainer(MGInfo<dim> &mg_info);

  /**
   * @brief Reinitialize the matrix-free object(s).
   */
  template <unsigned int degree, unsigned int quad_dim>
  void
  reinit(const dealii::Mapping<dim>                   &mapping,
         const DofHandler<dim>                        &dof_container,
         const ConstraintHandler<dim, degree, number> &constraint_container,
         const dealii::Quadrature<quad_dim>           &quad);

  /**
   * @brief Getter function for the matrix-free object (shared ptr).
   */
  [[nodiscard]] std::shared_ptr<
    dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>
  get_matrix_free() const;

  /**
   * @brief Getter function for the multigrid matrix-free object (shared ptr).
   */
  [[nodiscard]] std::shared_ptr<
    dealii::MatrixFree<dim, float, dealii::VectorizedArray<float>>>
  get_mg_matrix_free(unsigned int level) const;

private:
  /**
   * @brief Matrix-free object handler for non-multigrid data.
   */
  MatrixFreeHandler<dim, number> matrix_free;

  /**
   * @brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<MatrixFreeHandler<dim, float>> multigrid_matrix_free;

  /**
   * @brief Min multigrid level.
   */
  unsigned int min_level = 0;

  /**
   * @brief Max multigrid level.
   */
  unsigned int max_level = 0;
};

PRISMS_PF_END_NAMESPACE
