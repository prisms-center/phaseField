// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/vectorization.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/constraint_manager.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/system_wide.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Typedef for dealii::MatrixFree.
 */
using dealii::MatrixFree;

/**
 * @brief Typedef for solution block vector.
 */
template <typename number>
using BlockVector = dealii::LinearAlgebra::distributed::BlockVector<number>;

/**
 * @brief Typedef for solution vector.
 */
template <typename number>
using SolutionVector = typename BlockVector<number>::BlockType;

/**
 * @brief Containers for matrix free objects
 */
template <unsigned int dim, typename number>
class MatrixFreeManager
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;

  /**
   * @brief Constructor.
   */
  MatrixFreeManager() = default;

  /**
   * @brief Reinit.
   * @pre dof_manager and constraint_manager are reinit
   */
  template <unsigned int degree>
  void
  reinit(const DoFManager<dim, degree>                &dof_manager,
         const ConstraintManager<dim, degree, number> &constraint_manager);

  [[nodiscard]] const MatrixFree<dim, number> &
  get_shared_matrix_free() const;

  [[nodiscard]] const MatrixFree<dim, number> &
  get_generic_matrix_free() const;

  [[nodiscard]] const std::vector<MatrixFree<dim, number>> &
  get_shared_matrix_free_levels() const;

  [[nodiscard]] const MatrixFree<dim, number> &
  get_mg_shared_matrix_free(unsigned int relative_level) const;

  [[nodiscard]] const std::vector<MatrixFree<dim, number>> &
  get_generic_matrix_free_levels() const;

  [[nodiscard]] const MatrixFree<dim, number> &
  get_mg_generic_matrix_free(unsigned int relative_level) const;

  [[nodiscard]] std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>>
  get_block_partitioners(const std::set<unsigned int> &field_indices) const;

  [[nodiscard]] std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>>
  get_mg_block_partitioners(const std::set<unsigned int> &field_indices,
                            unsigned int                  relative_level) const;

  void
  initialize_block_vector(BlockVector<number>          &vec,
                          const std::set<unsigned int> &field_indices) const;

  void
  initialize_mg_block_vector(BlockVector<number>          &vec,
                             const std::set<unsigned int> &field_indices,
                             unsigned int                  relative_level) const;

private:
  /**
   * @brief MatrixFree object for every field.
   */
  MatrixFree<dim, number> shared_matrix_free;

  /**
   * @brief Generic Matrix-free object with a scalar and vector entry.
   */
  MatrixFree<dim, number> generic_matrix_free;

  /**
   * @brief MatrixFree object on each level for every field.
   */
  std::vector<MatrixFree<dim, number>> shared_matrix_free_levels;

  /**
   * @brief Generic Matrix-free object with a scalar and vector entry on each level.
   */
  std::vector<MatrixFree<dim, number>> generic_matrix_free_levels;
};

template <unsigned int dim, typename number>
template <unsigned int degree>
void
MatrixFreeManager<dim, number>::reinit(
  const DoFManager<dim, degree>                &dof_manager,
  const ConstraintManager<dim, degree, number> &constraint_manager)
{
  using AdditionalData = typename MatrixFree<dim, number>::AdditionalData;
  // TODO: maybe make the flags determined by user dependencies. (Though I doubt this
  // is a significant source of runtime)
  static const AdditionalData additional_data(
    AdditionalData::TasksParallelScheme::partition_partition,
    0,
    dealii::update_values | dealii::update_gradients | dealii::update_hessians |
      dealii::update_JxW_values | dealii::update_quadrature_points);

  const std::array<dealii::DoFHandler<dim>, 2> &generic_dof_handlers =
    dof_manager.get_dof_handlers();

  {
    const std::array<dealii::AffineConstraints<number>, 2> &generic_constraints =
      constraint_manager.get_generic_constraints();

    // Reinit shared MatrixFree
    shared_matrix_free.reinit(SystemWide<dim, degree>::mapping,
                              dof_manager.get_field_dof_handlers(),
                              constraint_manager.get_field_constraints(),
                              dealii::QGaussLobatto<1>(degree + 1),
                              // should dim really be 1?
                              additional_data);

    // Reinit generic MatrixFree
    generic_matrix_free.reinit(SystemWide<dim, degree>::mapping,
                               std::vector<const dealii::DoFHandler<dim> *>(
                                 {&generic_dof_handlers[0], &generic_dof_handlers[1]}),
                               std::vector<const dealii::AffineConstraints<number> *>(
                                 {&generic_constraints[0], &generic_constraints[1]}),
                               dealii::QGaussLobatto<1>(degree +
                                                        1)); // should dim really be 1?
  }
  const unsigned int num_levels = dof_manager.has_mg() ? dof_manager.num_levels() : 0;
  shared_matrix_free_levels.resize(num_levels);
  generic_matrix_free_levels.resize(num_levels);
  for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
    {
      const unsigned int       fine_level = num_levels - 1;
      const unsigned int       level      = fine_level - relative_level;
      MatrixFree<dim, number> &shared_mg_matrix_free =
        shared_matrix_free_levels[relative_level];
      MatrixFree<dim, number> &generic_mg_matrix_free =
        generic_matrix_free_levels[relative_level];

      const std::array<dealii::AffineConstraints<number>, 2> &generic_constraints =
        constraint_manager.get_mg_generic_constraints(relative_level);

      AdditionalData shared_additional_data = additional_data;
      shared_additional_data.mg_level       = level;

      // Reinit shared MatrixFree
      shared_mg_matrix_free.reinit(
        SystemWide<dim, degree>::mapping,
        dof_manager.get_field_dof_handlers(),
        constraint_manager.get_mg_field_constraints(relative_level),
        dealii::QGaussLobatto<1>(degree + 1), // should dim really be 1?
        shared_additional_data);

      AdditionalData generic_additional_data;
      generic_additional_data.mg_level = level;

      // Reinit generic MatrixFree
      generic_mg_matrix_free.reinit(
        SystemWide<dim, degree>::mapping,
        std::vector<const dealii::DoFHandler<dim> *>(
          {&generic_dof_handlers[0], &generic_dof_handlers[1]}),
        std::vector<const dealii::AffineConstraints<number> *>(
          {&generic_constraints[0], &generic_constraints[1]}),
        dealii::QGaussLobatto<1>(degree + 1), // should dim really be 1?
        generic_additional_data);
    }
}

template <unsigned int dim, typename number>
const MatrixFree<dim, number> &
MatrixFreeManager<dim, number>::get_shared_matrix_free() const
{
  return shared_matrix_free;
}

template <unsigned int dim, typename number>
const MatrixFree<dim, number> &
MatrixFreeManager<dim, number>::get_generic_matrix_free() const
{
  return generic_matrix_free;
}

template <unsigned int dim, typename number>
const std::vector<MatrixFree<dim, number>> &
MatrixFreeManager<dim, number>::get_shared_matrix_free_levels() const
{
  return shared_matrix_free_levels;
}

template <unsigned int dim, typename number>
const MatrixFree<dim, number> &
MatrixFreeManager<dim, number>::get_mg_shared_matrix_free(
  unsigned int relative_level) const
{
  return shared_matrix_free_levels[relative_level];
}

template <unsigned int dim, typename number>
const std::vector<MatrixFree<dim, number>> &
MatrixFreeManager<dim, number>::get_generic_matrix_free_levels() const
{
  return generic_matrix_free_levels;
}

template <unsigned int dim, typename number>
const MatrixFree<dim, number> &
MatrixFreeManager<dim, number>::get_mg_generic_matrix_free(
  unsigned int relative_level) const
{
  return generic_matrix_free_levels[relative_level];
}

template <unsigned int dim, typename number>
std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>>
MatrixFreeManager<dim, number>::get_block_partitioners(
  const std::set<unsigned int> &field_indices) const
{
  // These partitioners basically just provide the number of elements in a distributed way
  std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>> partitioners;
  partitioners.reserve(field_indices.size());
  for (unsigned int field_index : field_indices)
    {
      partitioners.push_back(shared_matrix_free.get_vector_partitioner(field_index));
    }
  return partitioners;
}

template <unsigned int dim, typename number>
std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>>
MatrixFreeManager<dim, number>::get_mg_block_partitioners(
  const std::set<unsigned int> &field_indices,
  unsigned int                  relative_level) const
{
  // These partitioners basically just provide the number of elements in a distributed way
  std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>> partitioners;
  partitioners.reserve(field_indices.size());
  for (unsigned int field_index : field_indices)
    {
      partitioners.push_back(
        shared_matrix_free_levels[relative_level].get_vector_partitioner(field_index));
    }
  return partitioners;
}

template <unsigned int dim, typename number>
void
MatrixFreeManager<dim, number>::initialize_block_vector(
  BlockVector<number>          &vec,
  const std::set<unsigned int> &field_indices) const
{
  vec.reinit(get_block_partitioners(field_indices));
}

template <unsigned int dim, typename number>
void
MatrixFreeManager<dim, number>::initialize_mg_block_vector(
  BlockVector<number>          &vec,
  const std::set<unsigned int> &field_indices,
  unsigned int                  relative_level) const
{
  vec.reinit(get_mg_block_partitioners(field_indices, relative_level));
}

PRISMS_PF_END_NAMESPACE
