// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/types.h>
#include <deal.II/base/vectorization.h>

#include <prismspf/core/exceptions.h>
#include <prismspf/core/field_container.h>

#include <prismspf/solvers/mf_operator.h>

#include "prismspf/core/group_solution_handler.h"

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::compute_operator(BlockVector<number>       &dst,
                                                  const BlockVector<number> &src) const
{
  data->cell_loop(&MFOperator::compute_local_operator, this, dst, src, true);
  if (scale_by_diagonal)
    {
      for (unsigned int block_index = 0; block_index < dst.n_blocks(); block_index++)
        {
          dst.block(block_index).scale(*(scaling_diagonal[block_index]));
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::compute_local_operator(
  const MatrixFree<dim, number>               &_data,
  BlockVector<number>                         &dst,
  const BlockVector<number>                   &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Construct FEEvaluation objects
  // The reason this is constructed here, rather than as a private member is because
  // compute_local_rhs is called by cell_loop, which multithreads. There would be data
  // races.
  FieldContainer<dim, degree, number> variable_list(field_attributes,
                                                    *solution_indexer,
                                                    relative_level,
                                                    dependency_map,
                                                    solve_block,
                                                    _data);

  // Initialize, evaluate, and submit based on user function.
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaluation flags for each variable
      variable_list.reinit_and_eval(cell, &src, read_plain);

      // Evaluate the user-defined pde at each quadrature point
      for (unsigned int quad = 0; quad < variable_list.get_n_q_points(); ++quad)
        {
          variable_list.set_q_point(quad);
          // Evaluate the function pointer (the user-defined pde)
          try
            {
              (pde_operator->*pde_op)(variable_list, *sim_timer, solve_block.id);
            }
          catch (...)
            {
              std::cerr << "Error: Exception thrown in equations during solve block "
                        << solve_block.id << "!" << std::endl;
              throw;
            }
        }

      // Integrate and add to global vector dst

      variable_list.integrate_and_distribute(&dst);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::compute_diagonal(BlockVector<number>       &dst,
                                                  const BlockVector<number> &src) const
{
  dst.reinit(src);
  data->cell_loop(&MFOperator::compute_local_diagonal, this, dst, src);

  // This is to make sure the preconditioner doesn't break down when there are
  // Dirichlet conditions (which lead to zero diagonal entries). The actual value
  // of these entries doesn't matter since they get overwritten by the constraints
  // application in the solver.
  set_zero_entries_to_one(dst);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::compute_local_diagonal(
  const MatrixFree<dim, number>               &_data,
  BlockVector<number>                         &diagonal,
  const BlockVector<number>                   &dummy_src, // just needs right shape
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  FieldContainer<dim, degree, number> variable_list(field_attributes,
                                                    *solution_indexer,
                                                    relative_level,
                                                    dependency_map,
                                                    solve_block,
                                                    _data);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Reinit the cell for all the dependencies
      variable_list.reinit_and_eval(cell, &dummy_src, false);

      // To get the diagonal of the "matrix", repeatedly "multiply the matrix" by a test x
      // vector (src vector) and store the solution in the i-th position in the y
      // vector (the diagonal).
      // First (i=1) x vector: I 0 0 0
      // Second(i=2) x vector: 0 I 0 0 ...
      // Submit zeros for everything except the diagonals
      for (unsigned int field_index : solve_block.field_indices)
        {
          if (field_attributes[field_index].field_type == TensorRank::Scalar)
            {
              compute_local_field_diagonal<TensorRank::Scalar>(variable_list,
                                                               diagonal,
                                                               field_index);
            }
          else if (field_attributes[field_index].field_type == TensorRank::Vector)
            {
              compute_local_field_diagonal<TensorRank::Vector>(variable_list,
                                                               diagonal,
                                                               field_index);
            }
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
void
MFOperator<dim, degree, number>::compute_local_field_diagonal(
  FieldContainer<dim, degree, number> &variable_list,
  BlockVector<number>                 &diagonal,
  unsigned int                         field_index) const
{
  // Number of nodes in the cell
  constexpr static unsigned int dofs_per_component =
    FieldContainer<dim, degree, number>::dofs_per_component;
  // "zero lhs vector"
  for (unsigned int some_field_index : solve_block.field_indices)
    {
      for (unsigned int i = 0; i < dofs_per_component; ++i)
        {
          variable_list.submit_dof_value(some_field_index, zero<Rank>(), i);
        }
    }
  // Object to hold the local diagonal
  dealii::AlignedVector<Value<Rank>> cell_diagonal(dofs_per_component, zero<Rank>());
  for (unsigned int i = 0; i < dofs_per_component; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_component; ++j)
        {
          variable_list.submit_dof_value(field_index,
                                         i == j ? identity<Rank>() : zero<Rank>(),
                                         j);
        }

      // Evaluate the dependencies based on the flags
      variable_list.eval_without_read();

      // Evaluate at each quadrature point
      for (unsigned int quad = 0; quad < variable_list.get_n_q_points(); ++quad)
        {
          variable_list.set_q_point(quad);
          try
            {
              (pde_operator->*pde_op)(variable_list, *sim_timer, solve_block.id);
            }
          catch (...)
            {
              std::cerr << "Error: Exception thrown in equations during solve block "
                        << solve_block.id << "!" << std::endl;
              throw;
            }
        }
      // Integrate the diagonal
      variable_list.integrate();

      // variable_list.eval(); //needed?
      variable_list.get_dof_value_to(cell_diagonal[i], field_index, i);
    }
  // Use FieldContainer interface to write local diagonal to global
  for (unsigned int i = 0; i < dofs_per_component; ++i)
    {
      variable_list.submit_dof_value(field_index, cell_diagonal[i], i);
    }
  variable_list.distribute(field_index, &diagonal);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::reinit_matrix_diagonal(const BlockVector<number> &shape)
{
  diagonal_entries->get_vector().reinit(shape);
  inverse_diagonal_entries->get_vector().reinit(shape);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::eval_matrix_diagonal()
{
  auto &diag     = diagonal_entries->get_vector();
  auto &inv_diag = inverse_diagonal_entries->get_vector();
  compute_diagonal(diag, diag);
  for (unsigned int i : diag.locally_owned_elements())
    {
      number diag_el = diag[i];
      inv_diag[i]    = 1.0 / diag_el;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::types::global_dof_index
MFOperator<dim, degree, number>::m() const
{
  Assert(data != nullptr, dealii::ExcNotInitialized());
  unsigned int sum = 0;
  for (unsigned int block_index = 0; block_index < solve_block.field_indices.size();
       ++block_index)
    {
      sum += data->get_vector_partitioner(block_index)->size();
    }
  return sum;
}

template <unsigned int dim, unsigned int degree, typename number>
number
MFOperator<dim, degree, number>::el([[maybe_unused]] const unsigned int &row,
                                    [[maybe_unused]] const unsigned int &col) const
{
  AssertThrow(false, FeatureNotImplemented("el()"));
  return 0.0;
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::clear()
{
  data = nullptr;
  diagonal_entries.reset();
  inverse_diagonal_entries.reset();
}

template <unsigned int dim, unsigned int degree, typename number>
const MatrixFree<dim, number> *
MFOperator<dim, degree, number>::get_matrix_free() const
{
  return data;
}

template <unsigned int dim, unsigned int degree, typename number>
const std::shared_ptr<dealii::DiagonalMatrix<BlockVector<number>>> &
MFOperator<dim, degree, number>::get_matrix_diagonal_inverse() const
{
  Assert(inverse_diagonal_entries.get() != nullptr && inverse_diagonal_entries->m() > 0,
         dealii::ExcNotInitialized());
  return inverse_diagonal_entries;
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::vmult(BlockVector<number>       &dst,
                                       const BlockVector<number> &src) const
{
  compute_operator(dst, src);
}

// NOLINTBEGIN(readability-identifier-naming)

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::Tvmult(BlockVector<number>       &dst,
                                        const BlockVector<number> &src) const
{
  vmult(dst, src);
}

// NOLINTEND(readability-identifier-naming)

#include "solvers/mf_operator.inst"

PRISMS_PF_END_NAMESPACE
