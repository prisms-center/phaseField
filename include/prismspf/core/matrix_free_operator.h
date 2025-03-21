// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/subscriptor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This is the abstract base class for the matrix-free implementation of some PDE.
 *
 * \tparam dim The number of dimensions in the problem.
 * \tparam degree The polynomial degree of the shape functions.
 * \tparam number Datatype to use for `LinearAlgebra::distributed::Vector<number>`. Either
 * double or float.
 */
template <int dim, int degree, typename number>
class matrixFreeOperator : public dealii::Subscriptor
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;
  using value_type = number;
  using size_type  = dealii::VectorizedArray<number>;

  /**
   * \brief Default constructor.
   */
  matrixFreeOperator() = default;

  /**
   * \brief Constructor for concurrent solves.
   */
  matrixFreeOperator(const userInputParameters<dim>                   &_user_inputs,
                     const std::map<unsigned int, variableAttributes> &_attributes_list);

  /**
   * \brief Constructor for single solves.
   */
  matrixFreeOperator(const userInputParameters<dim>                   &_user_inputs,
                     const unsigned int                               &_current_index,
                     const std::map<unsigned int, variableAttributes> &_attributes_list);

  /**
   * \brief Initialize operator.
   */
  void
  initialize(std::shared_ptr<const dealii::MatrixFree<dim, number, size_type>> _data,
             const std::vector<unsigned int> &selected_field_indexes =
               std::vector<unsigned int>());

  /**
   * \brief Return the number of DoFs.
   */
  dealii::types::global_dof_index
  m() const;

  /**
   * \brief Return the value of the matrix entry. This function is only valid when row ==
   * col and when the diagonal is initialized. Additionally, this is only used so that we
   * may compile. Trying to use this function will throw an error.
   */
  number
  el(const unsigned int &row, const unsigned int &col) const;

  /**
   * \brief Release all memory and return to state like having called the default
   * constructor.
   */
  void
  clear();

  /**
   * \brief Initialize a given vector with the MatrixFree object that this object
   * contains.
   */
  void
  initialize_dof_vector(VectorType &dst, unsigned int dof_handler_index = 0) const;

  /**
   * \brief Set constrained entries to one.
   */
  void
  set_constrained_entries_to_one(VectorType &dst) const;

  /**
   * \brief Get read access to the MatrixFree object stored with this operator.
   */
  std::shared_ptr<const dealii::MatrixFree<dim, number, size_type>>
  get_matrix_free() const;

  /**
   * \brief Get read access to the inverse diagonal of this operator.
   */
  const std::shared_ptr<dealii::DiagonalMatrix<VectorType>> &
  get_matrix_diagonal_inverse() const;

  /**
   * \brief Add the mappings from global to local solution vectors.
   */
  void
  add_global_to_local_mapping(
    const std::unordered_map<std::pair<unsigned int, dependencyType>,
                             unsigned int,
                             pairHash> &_global_to_local_solution);

  /**
   * \brief Add the solution subset for src vector.
   */
  void
  add_src_solution_subset(
    const std::vector<VectorType *> &_src_solution_subset = std::vector<VectorType *>());

  /**
   * \brief Matrix-vector multiplication.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * \brief Transpose matrix-vector multiplication.
   */
  void
  Tvmult(VectorType &dst, const VectorType &src) const;

  /**
   * \brief Compute the explicit update.
   */
  void
  compute_explicit_update(std::vector<VectorType *>       &dst,
                          const std::vector<VectorType *> &src) const;

  /**
   * \brief Compute the explicit update for postprocessed fields.
   */
  void
  compute_postprocess_explicit_update(std::vector<VectorType *>       &dst,
                                      const std::vector<VectorType *> &src) const;

  /**
   * \brief Compute a nonexplicit auxiliary update.
   */
  void
  compute_nonexplicit_auxiliary_update(std::vector<VectorType *>       &dst,
                                       const std::vector<VectorType *> &src) const;

  /**
   * \brief Compute the residual of this operator. This is the b in Ax=b.
   */
  void
  compute_residual(VectorType &dst, const VectorType &src) const;

  /**
   * \brief Compute the diagonal of this operator.
   */
  void
  compute_diagonal(unsigned int field_index);

  /**
   * \brief Get the user inputs (constant reference).
   */
  const userInputParameters<dim> &
  get_user_inputs() const;

  /**
   * \brief Get the index of the field that is currently being solved (copy).
   */
  unsigned int
  get_current_index() const;

  /**
   * \brief Get the timestep (copy).
   */
  number
  get_timestep() const;

protected:
  /**
   * \brief User-implemented class for the RHS of explicit equations.
   */
  virtual void
  compute_explicit_RHS(variableContainer<dim, degree, number> &variable_list,
                       const dealii::Point<dim, size_type>    &q_point_loc) const = 0;

  /**
   * \brief User-implemented class for the RHS of nonexplicit equations.
   */
  virtual void
  compute_nonexplicit_RHS(variableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, size_type>    &q_point_loc) const = 0;

  /**
   * \brief User-implemented class for the LHS of nonexplicit equations.
   */
  virtual void
  compute_nonexplicit_LHS(variableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, size_type>    &q_point_loc) const = 0;

  /**
   * \brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  virtual void
  compute_postprocess_explicit_RHS(
    variableContainer<dim, degree, number> &variable_list,
    const dealii::Point<dim, size_type>    &q_point_loc) const = 0;

private:
  /**
   * \brief Local computation of the explicit update.
   */
  void
  compute_local_explicit_update(
    const dealii::MatrixFree<dim, number, size_type> &data,
    std::vector<VectorType *>                        &dst,
    const std::vector<VectorType *>                  &src,
    const std::pair<unsigned int, unsigned int>      &cell_range) const;

  /**
   * \brief Local computation of the explicit update of postprocessed fields.
   */
  void
  compute_local_postprocess_explicit_update(
    const dealii::MatrixFree<dim, number, size_type> &data,
    std::vector<VectorType *>                        &dst,
    const std::vector<VectorType *>                  &src,
    const std::pair<unsigned int, unsigned int>      &cell_range) const;

  /**
   * \brief Local computation of the nonexplicit auxiliary update.
   */
  void
  compute_local_nonexplicit_auxiliary_update(
    const dealii::MatrixFree<dim, number, size_type> &data,
    std::vector<VectorType *>                        &dst,
    const std::vector<VectorType *>                  &src,
    const std::pair<unsigned int, unsigned int>      &cell_range) const;

  /**
   * \brief Local computation of the residual of the operator.
   */
  void
  compute_local_residual(const dealii::MatrixFree<dim, number, size_type> &data,
                         VectorType                                       &dst,
                         const VectorType                                 &src,
                         const std::pair<unsigned int, unsigned int> &cell_range) const;

  /**
   * \brief Local computation of the newton update of the operator.
   */
  void
  compute_local_newton_update(
    const dealii::MatrixFree<dim, number, size_type> &data,
    VectorType                                       &dst,
    const VectorType                                 &src,
    const std::pair<unsigned int, unsigned int>      &cell_range) const;

  /**
   * \brief Local computation of the diagonal of the operator.
   */
  void
  local_compute_diagonal(const dealii::MatrixFree<dim, number, size_type> &data,
                         VectorType                                       &dst,
                         const unsigned int                               &dummy,
                         const std::pair<unsigned int, unsigned int> &cell_range) const;

  /**
   * \brief The user-inputs.
   */
  const userInputParameters<dim> *user_inputs;

  /**
   * \brief The current index that is being solved.
   */
  unsigned int current_index = numbers::invalid_index;

  /**
   * \brief The attribute list of the relevant variables.
   */
  const std::map<unsigned int, variableAttributes> *attributes_list;

  /**
   * \brief Matrix-free object.
   */
  std::shared_ptr<const dealii::MatrixFree<dim, number, size_type>> data;

  /**
   * \brief Selected fields for which we'll evaluate.
   */
  std::vector<unsigned int> selected_fields;

  /**
   * \brief Indices of DoFs on edge in case the operator is used in GMG context.
   */
  std::vector<std::vector<unsigned int>> edge_constrained_indices;

  /**
   * \brief Mapping from global solution vectors to the local ones
   */
  std::unordered_map<std::pair<unsigned int, dependencyType>, unsigned int, pairHash>
    global_to_local_solution;

  /**
   * \brief Subset of fields that are necessary for the source.
   */
  std::vector<VectorType *> src_solution_subset;

  /**
   * \brief The diagonal matrix.
   */
  std::shared_ptr<dealii::DiagonalMatrix<VectorType>> diagonal_entries;

  /**
   * \brief The inverse diagonal matrix.
   */
  std::shared_ptr<dealii::DiagonalMatrix<VectorType>> inverse_diagonal_entries;
};

template <int dim, int degree, typename number>
matrixFreeOperator<dim, degree, number>::matrixFreeOperator(
  const userInputParameters<dim>                   &_user_inputs,
  const std::map<unsigned int, variableAttributes> &_attributes_list)
  : Subscriptor()
  , user_inputs(&_user_inputs)
  , attributes_list(&_attributes_list)
{}

template <int dim, int degree, typename number>
matrixFreeOperator<dim, degree, number>::matrixFreeOperator(
  const userInputParameters<dim>                   &_user_inputs,
  const unsigned int                               &_current_index,
  const std::map<unsigned int, variableAttributes> &_attributes_list)
  : Subscriptor()
  , user_inputs(&_user_inputs)
  , current_index(_current_index)
  , attributes_list(&_attributes_list)
{}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::initialize(
  std::shared_ptr<const dealii::MatrixFree<dim, number, size_type>> _data,
  const std::vector<unsigned int> &selected_field_indexes)
{
  data = _data;

  selected_fields.clear();
  if (selected_field_indexes.empty())
    {
      for (unsigned int i = 0; i < _data->n_components(); ++i)
        {
          selected_fields.push_back(i);
        }
    }
  else
    {
      for (unsigned int i = 0; i < selected_field_indexes.size(); ++i)
        {
          AssertIndexRange(selected_field_indexes[i], _data->n_components());
          for (unsigned int j = 0; j < selected_field_indexes.size(); ++j)
            {
              if (j != i)
                {
                  Assert(selected_field_indexes[j] != selected_field_indexes[i],
                         dealii::ExcMessage("Given row indices must be unique"));
                }
            }
          selected_fields.push_back(selected_field_indexes[i]);
        }
    }

  edge_constrained_indices.clear();
  edge_constrained_indices.resize(selected_fields.size());
}

template <int dim, int degree, typename number>
dealii::types::global_dof_index
matrixFreeOperator<dim, degree, number>::m() const
{
  Assert(data.get() != nullptr, dealii::ExcNotInitialized());
  unsigned int total_size = 0;
  for (const unsigned int field : selected_fields)
    {
      total_size += data->get_vector_partitioner(field)->size();
    }
  return total_size;
}

template <int dim, int degree, typename number>
number
matrixFreeOperator<dim, degree, number>::el(
  [[maybe_unused]] const unsigned int &row,
  [[maybe_unused]] const unsigned int &col) const
{
  AssertThrow(false, FeatureNotImplemented("el()"));
  return 0.0;
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::clear()
{
  data.reset();
  inverse_diagonal_entries.reset();
  global_to_local_solution.clear();
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::initialize_dof_vector(
  VectorType  &dst,
  unsigned int dof_handler_index) const
{
  data->initialize_dof_vector(dst, dof_handler_index);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::set_constrained_entries_to_one(
  VectorType &dst) const
{
  for (unsigned int j = 0; j < dealii::MatrixFreeOperators::BlockHelper::n_blocks(dst);
       ++j)
    {
      const std::vector<unsigned int> &constrained_dofs =
        data->get_constrained_dofs(selected_fields[j]);
      for (const auto constrained_dof : constrained_dofs)
        {
          dealii::MatrixFreeOperators::BlockHelper::subblock(dst, j).local_element(
            constrained_dof) = 1.0;
        }
      for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
        {
          dealii::MatrixFreeOperators::BlockHelper::subblock(dst, j).local_element(
            edge_constrained_indices[j][i]) = 1.0;
        }
    }
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::add_global_to_local_mapping(
  const std::unordered_map<std::pair<unsigned int, dependencyType>,
                           unsigned int,
                           pairHash> &_global_to_local_solution)
{
  global_to_local_solution = _global_to_local_solution;
}

template <int dim, int degree, typename number>
std::shared_ptr<const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>
matrixFreeOperator<dim, degree, number>::get_matrix_free() const
{
  return data;
}

template <int dim, int degree, typename number>
const std::shared_ptr<
  dealii::DiagonalMatrix<typename matrixFreeOperator<dim, degree, number>::VectorType>> &
matrixFreeOperator<dim, degree, number>::get_matrix_diagonal_inverse() const
{
  Assert(inverse_diagonal_entries.get() != nullptr && inverse_diagonal_entries->m() > 0,
         dealii::ExcNotInitialized());
  return inverse_diagonal_entries;
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::add_src_solution_subset(
  const std::vector<VectorType *> &_src_solution_subset)
{
  src_solution_subset = _src_solution_subset;
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_explicit_update(
  std::vector<VectorType *>       &dst,
  const std::vector<VectorType *> &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(!dst.empty(), dealii::ExcMessage("The dst vector must not be empty"));
  Assert(!src.empty(), dealii::ExcMessage("The src vector must not be empty"));

  this->data->cell_loop(&matrixFreeOperator::compute_local_explicit_update,
                        this,
                        dst,
                        src,
                        true);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_postprocess_explicit_update(
  std::vector<VectorType *>       &dst,
  const std::vector<VectorType *> &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(!dst.empty(), dealii::ExcMessage("The dst vector must not be empty"));
  Assert(!src.empty(), dealii::ExcMessage("The src vector must not be empty"));

  this->data->cell_loop(&matrixFreeOperator::compute_local_postprocess_explicit_update,
                        this,
                        dst,
                        src,
                        true);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_nonexplicit_auxiliary_update(
  std::vector<VectorType *>       &dst,
  const std::vector<VectorType *> &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(!dst.empty(), dealii::ExcMessage("The dst vector must not be empty"));
  Assert(!src.empty(), dealii::ExcMessage("The src vector must not be empty"));

  this->data->cell_loop(&matrixFreeOperator::compute_local_nonexplicit_auxiliary_update,
                        this,
                        dst,
                        src,
                        true);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_residual(VectorType       &dst,
                                                          const VectorType &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(!src_solution_subset.empty(),
         dealii::ExcMessage("The src_solution_subset vector must not be empty"));
  Assert(dst.size() != 0,
         dealii::ExcMessage("The dst vector should not have size equal to 0"));
  Assert(src.size() != 0,
         dealii::ExcMessage("The src vector should not have size equal to 0"));

  this->data->cell_loop(&matrixFreeOperator::compute_local_residual,
                        this,
                        dst,
                        src,
                        true);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::vmult(VectorType       &dst,
                                               const VectorType &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(dst.size() != 0,
         dealii::ExcMessage("The dst vector should not have size equal to 0"));
  Assert(src.size() != 0,
         dealii::ExcMessage("The src vector should not have size equal to 0"));

  this->data->cell_loop(&matrixFreeOperator::compute_local_newton_update,
                        this,
                        dst,
                        src,
                        true);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::Tvmult(VectorType       &dst,
                                                const VectorType &src) const
{
  this->vmult(dst, src);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_local_explicit_update(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  std::vector<VectorType *>                                              &dst,
  const std::vector<VectorType *>                                        &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  variableContainer<dim, degree, number> variable_list(data,
                                                       *attributes_list,
                                                       global_to_local_solution,
                                                       solveType::EXPLICIT_RHS);

  // Initialize, evaluate, and submit based on user function.
  variable_list.eval_local_operator(
    [this](variableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, size_type>    &q_point_loc)
    {
      this->compute_explicit_RHS(var_list, q_point_loc);
    },
    dst,
    src,
    cell_range);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_local_postprocess_explicit_update(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  std::vector<VectorType *>                                              &dst,
  const std::vector<VectorType *>                                        &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  variableContainer<dim, degree, number> variable_list(data,
                                                       *attributes_list,
                                                       global_to_local_solution,
                                                       solveType::POSTPROCESS);

  // Initialize, evaluate, and submit based on user function.
  variable_list.eval_local_operator(
    [this](variableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, size_type>    &q_point_loc)
    {
      this->compute_postprocess_explicit_RHS(var_list, q_point_loc);
    },
    dst,
    src,
    cell_range);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_local_nonexplicit_auxiliary_update(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  std::vector<VectorType *>                                              &dst,
  const std::vector<VectorType *>                                        &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  variableContainer<dim, degree, number> variable_list(data,
                                                       *attributes_list,
                                                       global_to_local_solution,
                                                       solveType::NONEXPLICIT_RHS);

  // Initialize, evaluate, and submit based on user function.
  variable_list.eval_local_operator(
    [this](variableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, size_type>    &q_point_loc)
    {
      this->compute_nonexplicit_RHS(var_list, q_point_loc);
    },
    dst,
    src,
    cell_range);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_local_residual(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  VectorType                                                             &dst,
  [[maybe_unused]] const VectorType                                      &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  variableContainer<dim, degree, number> variable_list(data,
                                                       *attributes_list,
                                                       global_to_local_solution,
                                                       solveType::NONEXPLICIT_RHS);

  // Initialize, evaluate, and submit based on user function.
  variable_list.eval_local_operator(
    [this](variableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, size_type>    &q_point_loc)
    {
      this->compute_nonexplicit_RHS(var_list, q_point_loc);
    },
    dst,
    src_solution_subset,
    cell_range);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_local_newton_update(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  VectorType                                                             &dst,
  const VectorType                                                       &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  variableContainer<dim, degree, number> variable_list(data,
                                                       *attributes_list,
                                                       global_to_local_solution,
                                                       solveType::NONEXPLICIT_LHS);

  // Initialize, evaluate, and submit based on user function. Note that the src solution
  // subset must not include the src vector.
  variable_list.eval_local_operator(
    [this](variableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, size_type>    &q_point_loc)
    {
      this->compute_nonexplicit_LHS(var_list, q_point_loc);
    },
    dst,
    src,
    src_solution_subset,
    cell_range);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_diagonal(unsigned int field_index)
{
  inverse_diagonal_entries.reset(new dealii::DiagonalMatrix<VectorType>());
  VectorType &inverse_diagonal = inverse_diagonal_entries->get_vector();
  data->initialize_dof_vector(inverse_diagonal, field_index);
  unsigned int dummy = 0;
  data->cell_loop(&matrixFreeOperator::local_compute_diagonal,
                  this,
                  inverse_diagonal,
                  dummy);

  set_constrained_entries_to_one(inverse_diagonal);

  for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
    {
      Assert(inverse_diagonal.local_element(i) > 0.0,
             dealii::ExcMessage(
               "No diagonal entry in a positive definite operator should be zero"));
      inverse_diagonal.local_element(i) = 1.0 / inverse_diagonal.local_element(i);
    }
}

template <int dim, int degree, typename number>
const userInputParameters<dim> &
matrixFreeOperator<dim, degree, number>::get_user_inputs() const
{
  Assert(user_inputs != nullptr, dealii::ExcNotInitialized());
  return *user_inputs;
}

template <int dim, int degree, typename number>
unsigned int
matrixFreeOperator<dim, degree, number>::get_current_index() const
{
  Assert(current_index != numbers::invalid_index, dealii::ExcNotInitialized());
  return current_index;
}

template <int dim, int degree, typename number>
number
matrixFreeOperator<dim, degree, number>::get_timestep() const
{
  Assert(user_inputs != nullptr, dealii::ExcNotInitialized());
  return user_inputs->temporal_discretization.dt;
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::local_compute_diagonal(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  VectorType                                                             &dst,
  [[maybe_unused]] const unsigned int                                    &dummy,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  variableContainer<dim, degree, number> variable_list(data,
                                                       *attributes_list,
                                                       global_to_local_solution,
                                                       solveType::NONEXPLICIT_LHS);

  // Initialize, evaluate, and submit diagonal based on user function.
  variable_list.eval_local_diagonal(
    [this](variableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, size_type>    &q_point_loc)
    {
      this->compute_nonexplicit_LHS(var_list, q_point_loc);
    },
    dst,
    src_solution_subset,
    cell_range);
}

PRISMS_PF_END_NAMESPACE