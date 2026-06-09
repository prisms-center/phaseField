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
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/system_wide.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

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

  [[nodiscard]] const std::vector<MatrixFree<dim, number>> &
  get_shared_matrix_free_levels() const;

  [[nodiscard]] const MatrixFree<dim, number> &
  get_shared_matrix_free(unsigned int relative_level) const;

  [[nodiscard]] const std::vector<MatrixFree<dim, number>> &
  get_generic_matrix_free_levels() const;

  [[nodiscard]] const MatrixFree<dim, number> &
  get_generic_matrix_free(unsigned int relative_level) const;

  [[nodiscard]] std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>>
  get_block_partitioners(const std::set<unsigned int> &field_indices,
                         unsigned int                  relative_level = 0) const;

  void
  initialize_block_vector(BlockVector<number>          &vec,
                          const std::set<unsigned int> &field_indices,
                          unsigned int                  relative_level = 0) const;

private:
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

  const std::vector<std::array<dealii::DoFHandler<dim>, 2>> &generic_dof_handlers_levels =
    dof_manager.get_dof_handlers_levels();
  const std::vector<std::array<dealii::AffineConstraints<number>, 2>>
    &generic_constraints_levels = constraint_manager.get_generic_constraints_levels();

  const unsigned int num_levels = generic_dof_handlers_levels.size();
  shared_matrix_free_levels.resize(num_levels);
  generic_matrix_free_levels.resize(num_levels);
  for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
    {
      MatrixFree<dim, number> &shared_matrix_free =
        shared_matrix_free_levels[relative_level];
      MatrixFree<dim, number> &generic_matrix_free =
        generic_matrix_free_levels[relative_level];

      // Reinit shared MatrixFree
      shared_matrix_free.reinit(SystemWide<dim, degree>::mapping,
                                dof_manager.get_field_dof_handlers(relative_level),
                                constraint_manager.get_field_constraints(relative_level),
                                dealii::QGaussLobatto<1>(degree + 1),
                                // should dim really be 1?
                                additional_data);

      // Reinit generic MatrixFree
      generic_matrix_free.reinit(SystemWide<dim, degree>::mapping,
                                 std::vector<const dealii::DoFHandler<dim> *>(
                                   {&generic_dof_handlers_levels[relative_level][0],
                                    &generic_dof_handlers_levels[relative_level][1]}),
                                 std::vector<const dealii::AffineConstraints<number> *>(
                                   {&generic_constraints_levels[relative_level][0],
                                    &generic_constraints_levels[relative_level][1]}),
                                 SystemWide<dim, degree>::quadrature);
    }
}

template <unsigned int dim, typename number>
const std::vector<MatrixFree<dim, number>> &
MatrixFreeManager<dim, number>::get_shared_matrix_free_levels() const
{
  return shared_matrix_free_levels;
}

template <unsigned int dim, typename number>
const MatrixFree<dim, number> &
MatrixFreeManager<dim, number>::get_shared_matrix_free(unsigned int relative_level) const
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
MatrixFreeManager<dim, number>::get_generic_matrix_free(unsigned int relative_level) const
{
  return generic_matrix_free_levels[relative_level];
}

template <unsigned int dim, typename number>
std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>>
MatrixFreeManager<dim, number>::get_block_partitioners(
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
  const std::set<unsigned int> &field_indices,
  unsigned int                  relative_level) const
{
  vec.reinit(get_block_partitioners(field_indices, relative_level));
}

PRISMS_PF_END_NAMESPACE
