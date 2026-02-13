// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

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
  using SolutionVector = GroupSolutionHandler<dim, number>::SolutionVector;
  using MatrixFree     = GroupSolutionHandler<dim, number>::MatrixFree;
  using ScalarValue    = dealii::VectorizedArray<number>;
  using VectorValue    = dealii::Tensor<1, dim, ScalarValue>;

  explicit InvMManager(const DofManager<dim> &dof_manager,
                       bool                   _calculate_scalar,
                       bool                   _calculate_vector)
    : calculate_scalar(_calculate_scalar)
    , calculate_vector(_calculate_vector)
  {
    auto dof_handlers = dof_manager.get_dof_handlers();
    data.resize(dof_handlers.size());
    for (unsigned int i = 0; i < data.size(); ++i)
      {
        for (unsigned int rank = 0; rank < 2; ++rank)
          {
            data[i][rank].reinit(SystemWide<dim, degree>::mapping,
                                 dof_handlers[i][rank],
                                 dealii::AffineConstraints<number>(), // make member?
                                 SystemWide<dim, degree>::quadrature);
          }
      }
  }

  /**
   * @brief Initialize.
   */
  void
  initialize()
  {
    for (unsigned int i = 0; i < data.size(); ++i)
      {
        if (calculate_scalar)
          {
            data[i][0]->initialize_dof_vector(jxw_scalar);
            data[i][0]->initialize_dof_vector(invm_scalar);
          }
        if (calculate_vector)
          {
            data[i][1]->initialize_dof_vector(jxw_vector);
            data[i][1]->initialize_dof_vector(invm_vector);
          }
      }
  }

  /**
   * @brief Compute element volume for the triangulation
   */
  void
  compute_scalar_invm()
  {
    data->cell_loop(&InvMManager::compute_local_scalar, this, jxw_scalar, 0);
    invert(invm_scalar, jxw_scalar);
  }

  void
  compute_vector_invm()
  {
    data->cell_loop(&InvMManager::compute_local_vector, this, jxw_vector, 0);
    invert(invm_vector, jxw_vector);
  }

  void
  compute_local_scalar(const MatrixFree                            &_data,
                       SolutionVector                              &dst,
                       [[maybe_unused]] const int                  &src,
                       const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    dealii::FEEvaluation<dim, degree, degree + 1, 1, number, ScalarValue> fe_eval(_data);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval->first.evaluate(dealii::update_JxW_values);
        for (unsigned int quad = 0; quad < SystemWide<dim, degree>::quadrature.size();
             ++quad)
          {
            fe_eval.submit_value(fe_eval.JxW(quad));
          }
        fe_eval.integrate_scatter(dealii::EvaluationFlags::values, dst);
      }
  }

  void
  compute_local_vector(const MatrixFree                            &_data,
                       SolutionVector                              &dst,
                       [[maybe_unused]] const int                  &src,
                       const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    dealii::FEEvaluation<dim, degree, degree + 1, 1, number, VectorValue> fe_eval(_data);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval->first.evaluate(dealii::update_JxW_values);
        for (unsigned int quad = 0; quad < SystemWide<dim, degree>::quadrature.size();
             ++quad)
          {
            fe_eval.submit_value(fe_eval.JxW(quad) * one);
          }
        fe_eval.integrate_scatter(dealii::EvaluationFlags::values, dst);
      }
  }

  void
  compute_invm()
  {
    if (calculate_scalar)
      {
        compute_scalar_invm();
      }
    if (calculate_vector)
      {
        compute_vector_invm();
      }
  }

private:
  void
  invert(SolutionVector &dst, const SolutionVector &src) const
  {
    Assert(dst.size() == src.size(), dealii::ExcInternalError());
    for (unsigned int i = 0; i < src.locally_owned_size(); ++i)
      {
        dst[i] = 1.0 / src[i];
      }
  }

  /**
   * @brief Matrix-free object.
   */
  std::vector<std::array<MatrixFree, 2>> data;

  bool calculate_scalar = false;
  bool calculate_vector = false;

  /**
   * @brief Vector that stores element volumes
   */
  SolutionVector jxw_scalar;
  SolutionVector jxw_vector;
  SolutionVector invm_scalar;
  SolutionVector invm_vector;

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
