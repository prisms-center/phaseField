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
   * @param dof_manager The DoF manager to get the triangulation and DoF handlers from.
   * @param _calculate_scalar Whether to calculate the scalar invm (element volume).
   * @param _calculate_vector Whether to calculate the vector invm (for vector fields).
   */
  explicit InvMManager(const DoFManager<dim, degree>                &dof_manager,
                       const ConstraintManager<dim, degree, number> &constraint_manager,
                       bool                                          _calculate_scalar,
                       bool                                          _calculate_vector)
    : num_levels(dof_manager.get_dof_handlers().size())
    , calculate_scalar(_calculate_scalar)
    , calculate_vector(_calculate_vector)
  {
    data.resize(num_levels);
    jxw_scalar.resize(num_levels);
    invm_scalar.resize(num_levels);
    invm_sqrt_scalar.resize(num_levels);
    jxw_vector.resize(num_levels);
    invm_vector.resize(num_levels);
    invm_sqrt_vector.resize(num_levels);
    reinit(dof_manager, constraint_manager);
  }

  void
  reinit(const DoFManager<dim, degree>                &dof_manager,
         const ConstraintManager<dim, degree, number> &constraint_manager)
  {
    const std::vector<std::array<dealii::DoFHandler<dim>, 2>> &dof_handlers =
      dof_manager.get_dof_handlers();
    const std::vector<std::array<dealii::AffineConstraints<number>, 2>> &constraints =
      constraint_manager.get_generic_constraints();
    for (unsigned int i = 0; i < data.size(); ++i)
      {
        for (unsigned int rank = 0; rank < 2; ++rank)
          {
            data[i][rank].reinit(SystemWide<dim, degree>::mapping,
                                 dof_handlers[i][rank],
                                 constraints[i][rank],
                                 SystemWide<dim, degree>::quadrature);
          }
      }
  }

  /**
   * @brief Recompute the invm vectors.
   */
  void
  compute_invm()
  {
    initialize();
    if (calculate_scalar)
      {
        compute_scalar_invm();
      }
    if (calculate_vector)
      {
        compute_vector_invm();
      }
  }

  /**
   * @brief Get the invm vector for a given rank and level.
   *
   * @param rank The tensor rank of the field (scalar or vector).
   * @param relative_level The relative level to get the invm for.
   */
  const SolutionVector<number> &
  get_invm(TensorRank rank, unsigned int relative_level) const
  {
    Assert((rank == TensorRank::Scalar && calculate_scalar) ||
             (rank == TensorRank::Vector && calculate_vector),
           dealii::ExcInternalError("Requested invm that was not calculated"));
    if (rank == TensorRank::Scalar)
      {
        return invm_scalar[relative_level];
      }
    // else
    {
      return invm_vector[relative_level];
    }
  }

  /**
   * @brief Get the integrated jxw vector for a given rank and level.
   *
   * @param rank The tensor rank of the field (scalar or vector).
   * @param relative_level The relative level to get the jxw for.
   */
  const SolutionVector<number> &
  get_jxw(TensorRank rank, unsigned int relative_level) const
  {
    Assert((rank == TensorRank::Scalar && calculate_scalar) ||
             (rank == TensorRank::Vector && calculate_vector),
           dealii::ExcInternalError("Requested invm that was not calculated"));
    if (rank == TensorRank::Scalar)
      {
        return jxw_scalar[relative_level];
      }
    // else
    {
      return jxw_vector[relative_level];
    }
  }

  /**
   * @brief Get the integrated jxw vector for a given rank and level.
   *
   * @param rank The tensor rank of the field (scalar or vector).
   * @param relative_level The relative level to get the jxw for.
   */
  const SolutionVector<number> &
  get_invm_sqrt(TensorRank rank, unsigned int relative_level) const
  {
    Assert((rank == TensorRank::Scalar && calculate_scalar) ||
             (rank == TensorRank::Vector && calculate_vector),
           dealii::ExcInternalError("Requested invm that was not calculated"));
    if (rank == TensorRank::Scalar)
      {
        return invm_sqrt_scalar[relative_level];
      }
    // else
    {
      return invm_sqrt_vector[relative_level];
    }
  }

  std::vector<const SolutionVector<number> *>
  get_invm(const std::vector<FieldAttributes> &field_container,
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

  std::vector<const SolutionVector<number> *>
  get_jxw(const std::vector<FieldAttributes> &field_container,
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

  std::vector<const SolutionVector<number> *>
  get_invm_sqrt(const std::vector<FieldAttributes> &field_container,
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

private:
  /**
   * @brief Initialize.
   */
  void
  initialize()
  {
    for (unsigned int i = 0; i < num_levels; ++i)
      {
        if (calculate_scalar)
          {
            jxw_scalar[i].reinit(data[i][0].get_vector_partitioner());
            invm_scalar[i].reinit(data[i][0].get_vector_partitioner());
            invm_sqrt_scalar[i].reinit(data[i][0].get_vector_partitioner());
          }
        if (calculate_vector)
          {
            jxw_vector[i].reinit(data[i][1].get_vector_partitioner());
            invm_vector[i].reinit(data[i][1].get_vector_partitioner());
            invm_sqrt_vector[i].reinit(data[i][1].get_vector_partitioner());
          }
      }
  }

  /**
   * @brief Compute element volume for the triangulation
   */
  void
  compute_scalar_invm()
  {
    for (unsigned int i = 0; i < data.size(); ++i)
      {
        data[i][0].cell_loop(&InvMManager::compute_local_scalar, this, jxw_scalar[i], 0);
        invert(invm_scalar[i], jxw_scalar[i]);
        //
        sqrt(invm_sqrt_scalar[i], invm_scalar[i]);
        jxw_scalar[i].update_ghost_values();
        invm_scalar[i].update_ghost_values();
        invm_sqrt_scalar[i].update_ghost_values();
      }
  }

  void
  compute_vector_invm()
  {
    for (unsigned int i = 0; i < data.size(); ++i)
      {
        data[i][1].cell_loop(&InvMManager::compute_local_vector, this, jxw_vector[i], 0);
        invert(invm_vector[i], jxw_vector[i]);
        //
        sqrt(invm_sqrt_vector[i], invm_vector[i]);
        jxw_vector[i].update_ghost_values();
        invm_vector[i].update_ghost_values();
        invm_sqrt_vector[i].update_ghost_values();
      }
  }

  void
  compute_local_scalar(const MatrixFree<dim, number>               &_data,
                       SolutionVector<number>                      &dst,
                       [[maybe_unused]] const int                  &src,
                       const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    dealii::FEEvaluation<dim, degree, degree + 1, 1, number, ScalarValue> fe_eval(_data);
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
    dealii::FEEvaluation<dim, degree, degree + 1, dim, number> fe_eval(_data);
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

  /**
   * @brief Matrix-free object.
   */
  std::vector<std::array<MatrixFree<dim, number>, 2>> data;

  unsigned int num_levels;

  bool calculate_scalar = false;
  bool calculate_vector = false;

  /**
   * @brief Vector that stores element volumes
   */
  std::vector<SolutionVector<number>> jxw_scalar;
  std::vector<SolutionVector<number>> jxw_vector;
  std::vector<SolutionVector<number>> invm_scalar;
  std::vector<SolutionVector<number>> invm_vector;
  std::vector<SolutionVector<number>> invm_sqrt_scalar;
  std::vector<SolutionVector<number>> invm_sqrt_vector;

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

PRISMS_PF_END_NAMESPACE
