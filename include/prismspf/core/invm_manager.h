// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <prismspf/core/matrix_free_manager.h>
#include <prismspf/core/system_wide.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief A little class that computes the element volume for our triangulation.
 */
template <unsigned int dim, unsigned int degree, typename number>
class InvMManager
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;

  /**
   * @brief Constructor.
   */
  InvMManager() = default;

  /**
   * @brief Reserve space.
   */
  void
  reinit(const MatrixFreeManager<dim, number> &matrix_free_manager);

  /**
   * @brief Recompute the invm vectors.
   */
  void
  compute_invm();

  /**
   * @brief Get the invm vector for a given rank and level.
   *
   * @param rank The tensor rank of the field (scalar or vector).
   */
  const SolutionVector<number> &
  get_invm(TensorRank rank, unsigned int relative_level = -1) const;

  /**
   * @brief Get the integrated jxw vector for a given rank and level.
   *
   * @param rank The tensor rank of the field (scalar or vector).
   */
  const SolutionVector<number> &
  get_jxw(TensorRank rank, unsigned int relative_level = -1) const;

  /**
   * @brief Get the integrated jxw vector for a given rank.
   *
   * @param rank The tensor rank of the field (scalar or vector).
   */
  const SolutionVector<number> &
  get_invm_sqrt(TensorRank rank, unsigned int relative_level = -1) const;

  std::vector<const SolutionVector<number> *>
  get_invm(const std::vector<FieldAttributes> &field_container,
           const std::set<unsigned int>       &field_indices,
           unsigned int                        relative_level = -1) const;

  std::vector<const SolutionVector<number> *>
  get_jxw(const std::vector<FieldAttributes> &field_container,
          const std::set<unsigned int>       &field_indices,
          unsigned int                        relative_level = -1) const;

  std::vector<const SolutionVector<number> *>
  get_invm_sqrt(const std::vector<FieldAttributes> &field_container,
                const std::set<unsigned int>       &field_indices,
                unsigned int                        relative_level = -1) const;

private: // todo: move to outside class body
  const MatrixFreeManager<dim, number> &
  mf_manager()
  {
    Assert(mf_manager_ptr != nullptr, dealii::ExcNotInitialized());
    return *mf_manager_ptr;
  }

  /**
   * @brief Initialize.
   */
  void
  initialize()
  {
    const MatrixFree<dim, number> &matrix_free = mf_manager().get_generic_matrix_free();
    jxw_scalar.reinit(matrix_free.get_vector_partitioner(0));
    invm_scalar.reinit(matrix_free.get_vector_partitioner(0));
    invm_sqrt_scalar.reinit(matrix_free.get_vector_partitioner(0));

    jxw_vector.reinit(matrix_free.get_vector_partitioner(1));
    invm_vector.reinit(matrix_free.get_vector_partitioner(1));
    invm_sqrt_vector.reinit(matrix_free.get_vector_partitioner(1));
    for (unsigned int i = 0; i < num_levels; ++i)
      {
        const MatrixFree<dim, number> &mg_matrix_free =
          mf_manager().get_mg_generic_matrix_free(i);
        mg_jxw_scalar[i].reinit(mg_matrix_free.get_vector_partitioner(0));
        mg_invm_scalar[i].reinit(mg_matrix_free.get_vector_partitioner(0));
        mg_invm_sqrt_scalar[i].reinit(mg_matrix_free.get_vector_partitioner(0));

        mg_jxw_vector[i].reinit(mg_matrix_free.get_vector_partitioner(1));
        mg_invm_vector[i].reinit(mg_matrix_free.get_vector_partitioner(1));
        mg_invm_sqrt_vector[i].reinit(mg_matrix_free.get_vector_partitioner(1));
      }
  }

  /**
   * @brief Compute element volume for the triangulation
   */
  void
  compute_scalar_invm()
  {
    const MatrixFree<dim, number> &matrix_free = mf_manager().get_generic_matrix_free();
    matrix_free.cell_loop(&InvMManager::compute_local_scalar, this, jxw_scalar, 0);
    invert(invm_scalar, jxw_scalar);
    //
    sqrt(invm_sqrt_scalar, invm_scalar);
    jxw_scalar.update_ghost_values();
    invm_scalar.update_ghost_values();
    invm_sqrt_scalar.update_ghost_values();
    for (unsigned int i = 0; i < num_levels; ++i)
      {
        const MatrixFree<dim, number> &mg_matrix_free =
          mf_manager().get_mg_generic_matrix_free(i);
        mg_matrix_free.cell_loop(&InvMManager::compute_local_scalar,
                                 this,
                                 mg_jxw_scalar[i],
                                 0);
        invert(mg_invm_scalar[i], mg_jxw_scalar[i]);
        //
        sqrt(mg_invm_sqrt_scalar[i], mg_invm_scalar[i]);
        mg_jxw_scalar[i].update_ghost_values();
        mg_invm_scalar[i].update_ghost_values();
        mg_invm_sqrt_scalar[i].update_ghost_values();
      }
  }

  void
  compute_vector_invm()
  {
    const MatrixFree<dim, number> &matrix_free = mf_manager().get_generic_matrix_free();
    matrix_free.cell_loop(&InvMManager::compute_local_vector, this, jxw_vector, 0);
    invert(invm_vector, jxw_vector);
    //
    sqrt(invm_sqrt_vector, invm_vector);
    jxw_vector.update_ghost_values();
    invm_vector.update_ghost_values();
    invm_sqrt_vector.update_ghost_values();
    for (unsigned int i = 0; i < num_levels; ++i)
      {
        const MatrixFree<dim, number> &mg_matrix_free =
          mf_manager().get_mg_generic_matrix_free(i);
        mg_matrix_free.cell_loop(&InvMManager::compute_local_vector,
                                 this,
                                 mg_jxw_vector[i],
                                 0);
        invert(mg_invm_vector[i], mg_jxw_vector[i]);
        //
        sqrt(mg_invm_sqrt_vector[i], mg_invm_vector[i]);
        mg_jxw_vector[i].update_ghost_values();
        mg_invm_vector[i].update_ghost_values();
        mg_invm_sqrt_vector[i].update_ghost_values();
      }
  }

  void
  compute_local_scalar(const MatrixFree<dim, number>               &_data,
                       SolutionVector<number>                      &dst,
                       [[maybe_unused]] const int                  &src,
                       const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    dealii::FEEvaluation<dim, degree, degree + 1, 1, number> fe_eval(_data, 0);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        for (unsigned int quad = 0; quad < SystemWide<dim, degree>::quadrature.size();
             ++quad)
          {
            fe_eval.submit_value(dealii::make_vectorized_array<number>(1.0), quad);
          }
        fe_eval.integrate_scatter(dealii::EvaluationFlags::values, dst);
      }
  }

  void
  compute_local_vector(const MatrixFree<dim, number>               &_data,
                       SolutionVector<number>                      &dst,
                       [[maybe_unused]] const int                  &src,
                       const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    dealii::FEEvaluation<dim, degree, degree + 1, dim, number> fe_eval(_data, 1);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        for (unsigned int quad = 0; quad < SystemWide<dim, degree>::quadrature.size();
             ++quad)
          {
            fe_eval.submit_value(one, quad);
          }
        fe_eval.integrate_scatter(dealii::EvaluationFlags::values, dst);
      }
  }

  void
  invert(SolutionVector<number> &dst, const SolutionVector<number> &src) const
  {
    src.update_ghost_values();
    dst.update_ghost_values();
    Assert(dst.size() == src.size(), dealii::ExcInternalError());
    for (unsigned int i = 0; i < src.locally_owned_size(); ++i)
      {
        number src_el = src.local_element(i);
        // for some reason, some elements are zero
        dst.local_element(i) = src_el == 0 ? 0.0 : 1.0 / src_el;
      }
  }

  void
  sqrt(SolutionVector<number> &dst, const SolutionVector<number> &src) const
  {
    src.update_ghost_values();
    dst.update_ghost_values();
    Assert(dst.size() == src.size(), dealii::ExcInternalError());
    for (unsigned int i = 0; i < src.locally_owned_size(); ++i)
      {
        dst.local_element(i) = std::sqrt(src.local_element(i));
      }
  }

  unsigned int num_levels = 0;

  /**
   * @brief Vector that stores element volumes
   */
  SolutionVector<number>              jxw_scalar;
  SolutionVector<number>              jxw_vector;
  SolutionVector<number>              invm_scalar;
  SolutionVector<number>              invm_vector;
  SolutionVector<number>              invm_sqrt_scalar;
  SolutionVector<number>              invm_sqrt_vector;
  std::vector<SolutionVector<number>> mg_jxw_scalar;
  std::vector<SolutionVector<number>> mg_jxw_vector;
  std::vector<SolutionVector<number>> mg_invm_scalar;
  std::vector<SolutionVector<number>> mg_invm_vector;
  std::vector<SolutionVector<number>> mg_invm_sqrt_scalar;
  std::vector<SolutionVector<number>> mg_invm_sqrt_vector;

  const MatrixFreeManager<dim, number> *mf_manager_ptr = nullptr;

  inline static const VectorValue one = []()
  {
    VectorValue one1;
    for (unsigned int i = 0; i < dim; ++i)
      {
        one1[i] = 1.0;
      }
    return one1;
  }();
};

template <unsigned int dim, unsigned int degree, typename number>
inline void
InvMManager<dim, degree, number>::reinit(
  const MatrixFreeManager<dim, number> &matrix_free_manager)
{
  mf_manager_ptr = &matrix_free_manager;
  num_levels     = matrix_free_manager.get_generic_matrix_free_levels().size();
  mg_jxw_scalar.resize(num_levels);
  mg_invm_scalar.resize(num_levels);
  mg_invm_sqrt_scalar.resize(num_levels);
  mg_jxw_vector.resize(num_levels);
  mg_invm_vector.resize(num_levels);
  mg_invm_sqrt_vector.resize(num_levels);
  initialize();
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
InvMManager<dim, degree, number>::compute_invm()
{
  compute_scalar_invm();
  compute_vector_invm();
}

template <unsigned int dim, unsigned int degree, typename number>
inline const SolutionVector<number> &
InvMManager<dim, degree, number>::get_invm(TensorRank   rank,
                                           unsigned int relative_level) const
{
  if (relative_level == -1)
    {
      if (rank == TensorRank::Scalar)
        {
          return invm_scalar;
        }
      // else
      {
        return invm_vector;
      }
    }
  if (rank == TensorRank::Scalar)
    {
      return mg_invm_scalar[relative_level];
    }
  // else
  {
    return mg_invm_vector[relative_level];
  }
}

template <unsigned int dim, unsigned int degree, typename number>
inline const SolutionVector<number> &
InvMManager<dim, degree, number>::get_jxw(TensorRank   rank,
                                          unsigned int relative_level) const
{
  if (relative_level == -1)
    {
      if (rank == TensorRank::Scalar)
        {
          return jxw_scalar;
        }
      // else
      {
        return jxw_vector;
      }
    }
  if (rank == TensorRank::Scalar)
    {
      return mg_jxw_scalar[relative_level];
    }
  // else
  {
    return mg_jxw_vector[relative_level];
  }
}

template <unsigned int dim, unsigned int degree, typename number>
inline const SolutionVector<number> &
InvMManager<dim, degree, number>::get_invm_sqrt(TensorRank   rank,
                                                unsigned int relative_level) const
{
  if (relative_level == -1)
    {
      if (rank == TensorRank::Scalar)
        {
          return invm_sqrt_scalar;
        }
      // else
      {
        return invm_sqrt_vector;
      }
    }
  if (rank == TensorRank::Scalar)
    {
      return mg_invm_sqrt_scalar[relative_level];
    }
  // else
  {
    return mg_invm_sqrt_vector[relative_level];
  }
}

template <unsigned int dim, unsigned int degree, typename number>
inline std::vector<const SolutionVector<number> *>
InvMManager<dim, degree, number>::get_invm(
  const std::vector<FieldAttributes> &field_container,
  const std::set<unsigned int>       &field_indices,
  unsigned int                        relative_level) const
{
  const unsigned int                          num_blocks = field_indices.size();
  std::vector<const SolutionVector<number> *> out;
  out.reserve(num_blocks);
  for (unsigned int field_index : field_indices)
    {
      out.push_back(&get_invm(field_container[field_index].field_type, relative_level));
    }
  return out;
}

template <unsigned int dim, unsigned int degree, typename number>
inline std::vector<const SolutionVector<number> *>
InvMManager<dim, degree, number>::get_jxw(
  const std::vector<FieldAttributes> &field_container,
  const std::set<unsigned int>       &field_indices,
  unsigned int                        relative_level) const
{
  const unsigned int                          num_blocks = field_indices.size();
  std::vector<const SolutionVector<number> *> out;
  out.reserve(num_blocks);
  for (unsigned int field_index : field_indices)
    {
      out.push_back(&get_jxw(field_container[field_index].field_type, relative_level));
    }
  return out;
}

template <unsigned int dim, unsigned int degree, typename number>
inline std::vector<const SolutionVector<number> *>
InvMManager<dim, degree, number>::get_invm_sqrt(
  const std::vector<FieldAttributes> &field_container,
  const std::set<unsigned int>       &field_indices,
  unsigned int                        relative_level) const
{
  const unsigned int                          num_blocks = field_indices.size();
  std::vector<const SolutionVector<number> *> out;
  out.reserve(num_blocks);
  for (unsigned int field_index : field_indices)
    {
      out.push_back(
        &get_invm_sqrt(field_container[field_index].field_type, relative_level));
    }
  return out;
}

PRISMS_PF_END_NAMESPACE
