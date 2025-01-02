#include <deal.II/base/exceptions.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <core/solveTypeEnums.h>
#include <core/varTypeEnums.h>
#include <core/variableAttributes.h>
#include <core/variableContainer.h>
#include <string>

template <int dim, int degree, typename T>
variableContainer<dim, degree, T>::variableContainer(
  const dealii::MatrixFree<dim, double> &data,
  AttributesList                         _subset_attributes,
  solveType                              _solve_type)
  : subset_attributes(std::move(_subset_attributes))
  , solve_type(_solve_type)
{
  for (const auto &[var_index, variable] : subset_attributes)
    {
      uint      dof_no           = (solve_type != solveType::POSTPROCESS) ? var_index : 0;
      EvalFlags evaluation_flags = variable.eval_flags_for_solve(solve_type);
      bool      var_needed =
        (evaluation_flags != dealii::EvaluationFlags::nothing) || variable.is_pp;
      bool change_needed =
        (variable.eval_flags_change_nonexplicit_LHS != dealii::EvaluationFlags::nothing);
      if (var_needed)
        {
          if (variable.var_type == SCALAR)
            {
              scalar_vars_map.emplace(var_index,
                                      std::make_unique<scalar_FEEval>(data, dof_no));
            }
          else
            {
              vector_vars_map.emplace(var_index,
                                      std::make_unique<vector_FEEval>(data, dof_no));
            }
        }
      if (change_needed)
        {
          if (variable.var_type == SCALAR)
            {
              scalar_change_in_vars_map.emplace(var_index,
                                                std::make_unique<scalar_FEEval>(data,
                                                                                dof_no));
            }
          else
            {
              vector_change_in_vars_map.emplace(var_index,
                                                std::make_unique<vector_FEEval>(data,
                                                                                dof_no));
            }
        }
    }
}

template <int dim, int degree, typename T>
unsigned int
variableContainer<dim, degree, T>::get_num_q_points() const
{
  if (!scalar_vars_map.empty())
    {
      return scalar_vars_map.begin()->second->n_q_points;
    }
  else if (!vector_vars_map.empty())
    {
      return vector_vars_map.begin()->second->n_q_points;
    }
  else if (!scalar_change_in_vars_map.empty())
    {
      return scalar_change_in_vars_map.begin()->second->n_q_points;
    }
  else if (!vector_change_in_vars_map.empty())
    {
      return vector_change_in_vars_map.begin()->second->n_q_points;
    }
  else
    {
      AssertThrow(false,
                  dealii::ExcMessage(
                    "PRISMS-PF Error: When trying to access the number of quadrature "
                    "points, all FEEvaluation object containers were empty."));
    }
}

template <int dim, int degree, typename T>
dealii::Point<dim, T>
variableContainer<dim, degree, T>::get_q_point_location() const
{
  if (!scalar_vars_map.empty())
    {
      return scalar_vars_map.begin()->second->quadrature_point(q_point);
    }
  else if (!vector_vars_map.empty())
    {
      return vector_vars_map.begin()->second->quadrature_point(q_point);
    }
  else if (!scalar_change_in_vars_map.empty())
    {
      return scalar_change_in_vars_map.begin()->second->quadrature_point(q_point);
    }
  else if (!vector_change_in_vars_map.empty())
    {
      return vector_change_in_vars_map.begin()->second->quadrature_point(q_point);
    }
  else
    {
      AssertThrow(false,
                  dealii::ExcMessage(
                    "PRISMS-PF Error: When trying to access the quadrature point "
                    "location, all FEEvaluation object containers were empty."));
    }
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::reinit_and_eval(const std::vector<vectorType *> &src,
                                                   unsigned int                     cell)
{
  for (const auto &[var_index, variable] : subset_attributes)
    {
      if (variable.var_needed(solve_type) || variable.is_pp)
        {
          const uint dof_no = var_index;
          if (variable.var_type == SCALAR)
            {
              auto *scalar_FEEval_ptr = scalar_vars_map[var_index].get();
              scalar_FEEval_ptr->reinit(cell);
              scalar_FEEval_ptr->read_dof_values(*src[dof_no]);
              scalar_FEEval_ptr->evaluate(variable.eval_flags_for_solve(solve_type));
            }
          else
            {
              auto *vector_FEEval_ptr = vector_vars_map[var_index].get();
              vector_FEEval_ptr->reinit(cell);
              vector_FEEval_ptr->read_dof_values(*src[dof_no]);
              vector_FEEval_ptr->evaluate(variable.eval_flags_for_solve(solve_type));
            }
        }
    }
}

/**
 * This is specialized for the LHS where a change in solution is needed. The RHS method
 * takes the src as a vector of vectorTypes.
 */
template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::reinit_and_eval_change_in_solution(
  const vectorType &src,
  unsigned int      cell,
  const uint       &var_index)
{
  const variableAttributes variable = subset_attributes.at(var_index);
  if (variable.var_type == SCALAR)
    {
      auto *scalar_FEEval_ptr = scalar_change_in_vars_map[var_index].get();
      scalar_FEEval_ptr->reinit(cell);
      scalar_FEEval_ptr->read_dof_values(src);
      scalar_FEEval_ptr->evaluate(variable.eval_flags_residual_nonexplicit_LHS);
    }
  else
    {
      auto *vector_FEEval_ptr = vector_change_in_vars_map[var_index].get();
      vector_FEEval_ptr->reinit(cell);
      vector_FEEval_ptr->read_dof_values(src);
      vector_FEEval_ptr->evaluate(variable.eval_flags_residual_nonexplicit_LHS);
    }
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::reinit(unsigned int cell)
{
  for (const auto &[var_index, variable] : subset_attributes)
    {
      if (variable.var_needed(solve_type) || variable.is_pp)
        {
          if (variable.var_type == SCALAR)
            {
              scalar_vars_map[var_index]->reinit(cell);
            }
          else
            {
              vector_vars_map[var_index]->reinit(cell);
            }
        }
    }
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::integrate_and_distribute(
  std::vector<vectorType *> &dst)
{
  for (const auto &[var_index, variable] : subset_attributes)
    {
      if (variable.var_needed(solve_type) || variable.is_pp)
        {
          const uint dof_no = var_index;
          if (variable.var_type == SCALAR)
            {
              auto *scalar_FEEval_ptr = scalar_vars_map[var_index].get();
              scalar_FEEval_ptr->integrate(variable.residual_flags_for_solve(solve_type));
              scalar_FEEval_ptr->distribute_local_to_global(*dst[dof_no]);
            }
          else
            {
              auto *vector_FEEval_ptr = vector_vars_map[var_index].get();
              vector_FEEval_ptr->integrate(variable.residual_flags_for_solve(solve_type));
              vector_FEEval_ptr->distribute_local_to_global(*dst[dof_no]);
            }
        }
    }
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::integrate_and_distribute_change_in_solution_LHS(
  vectorType &dst,
  const uint &var_index)
{
  // integrate
  const variableAttributes variable = subset_attributes.at(var_index);
  if (variable.var_type == SCALAR)
    {
      auto *scalar_FEEval_ptr = scalar_change_in_vars_map[var_index].get();
      scalar_FEEval_ptr->integrate(variable.residual_flags_for_solve(solve_type));
      scalar_FEEval_ptr->distribute_local_to_global(dst);
    }
  else
    {
      auto *vector_FEEval_ptr = vector_change_in_vars_map[var_index].get();
      vector_FEEval_ptr->integrate(variable.residual_flags_for_solve(solve_type));
      vector_FEEval_ptr->distribute_local_to_global(dst);
    }
}

// Need to add index checking to these functions so that an error is thrown if
// the index wasn't set
template <int dim, int degree, typename T>
T
variableContainer<dim, degree, T>::get_scalar_value(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, SCALAR, dealii::EvaluationFlags::values, false);
  return scalar_vars_map.at(global_variable_index)->get_value(q_point);
}

template <int dim, int degree, typename T>
dealii::Tensor<1, dim, T>
variableContainer<dim, degree, T>::get_scalar_gradient(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, SCALAR, dealii::EvaluationFlags::gradients, false);
  return scalar_vars_map.at(global_variable_index)->get_gradient(q_point);
}

template <int dim, int degree, typename T>
dealii::Tensor<2, dim, T>
variableContainer<dim, degree, T>::get_scalar_hessian(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, SCALAR, dealii::EvaluationFlags::hessians, false);
  return scalar_vars_map.at(global_variable_index)->get_hessian(q_point);
}

template <int dim, int degree, typename T>
dealii::Tensor<1, dim, T>
variableContainer<dim, degree, T>::get_vector_value(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, VECTOR, dealii::EvaluationFlags::values, false);
  return vector_vars_map.at(global_variable_index)->get_value(q_point);
}

template <int dim, int degree, typename T>
dealii::Tensor<2, dim, T>
variableContainer<dim, degree, T>::get_vector_gradient(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, VECTOR, dealii::EvaluationFlags::gradients, false);
  return vector_vars_map.at(global_variable_index)->get_gradient(q_point);
}

template <int dim, int degree, typename T>
dealii::Tensor<3, dim, T>
variableContainer<dim, degree, T>::get_vector_hessian(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, VECTOR, dealii::EvaluationFlags::hessians, false);
  return vector_vars_map.at(global_variable_index)->get_hessian(q_point);
}

template <int dim, int degree, typename T>
T
variableContainer<dim, degree, T>::get_change_in_scalar_value(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, SCALAR, dealii::EvaluationFlags::values, true);
  return scalar_change_in_vars_map.at(global_variable_index)->get_value(q_point);
}

template <int dim, int degree, typename T>
dealii::Tensor<1, dim, T>
variableContainer<dim, degree, T>::get_change_in_scalar_gradient(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, SCALAR, dealii::EvaluationFlags::gradients, true);
  return scalar_change_in_vars_map.at(global_variable_index)->get_gradient(q_point);
}

template <int dim, int degree, typename T>
dealii::Tensor<2, dim, T>
variableContainer<dim, degree, T>::get_change_in_scalar_hessian(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, SCALAR, dealii::EvaluationFlags::hessians, true);
  return scalar_change_in_vars_map.at(global_variable_index)->get_hessian(q_point);
}

template <int dim, int degree, typename T>
dealii::Tensor<1, dim, T>
variableContainer<dim, degree, T>::get_change_in_vector_value(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, VECTOR, dealii::EvaluationFlags::values, true);
  return vector_change_in_vars_map.at(global_variable_index)->get_value(q_point);
}

template <int dim, int degree, typename T>
dealii::Tensor<2, dim, T>
variableContainer<dim, degree, T>::get_change_in_vector_gradient(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, VECTOR, dealii::EvaluationFlags::gradients, true);
  return vector_change_in_vars_map.at(global_variable_index)->get_gradient(q_point);
}

template <int dim, int degree, typename T>
dealii::Tensor<3, dim, T>
variableContainer<dim, degree, T>::get_change_in_vector_hessian(
  unsigned int global_variable_index) const
{
  AssertValid(global_variable_index, VECTOR, dealii::EvaluationFlags::hessians, true);
  return vector_change_in_vars_map.at(global_variable_index)->get_hessian(q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::AssertValid(const uint      &var_index,
                                               const fieldType &field_type,
                                               const EvalFlags &eval_flag,
                                               bool             is_change) const
{
#ifdef DEBUG
  std::string field_string = (field_type == SCALAR) ? "scalar" : "vector";
  Assert((subset_attributes.find(var_index) != subset_attributes.end()) &&
           bool((!is_change
                   ? subset_attributes.at(var_index).eval_flags_for_solve(solve_type)
                   : subset_attributes.at(var_index).eval_flags_change_nonexplicit_LHS) &
                eval_flag),
         dealii::ExcMessage(
           "PRISMS-PF Error: Attempted access of a variable value that was not marked as "
           "needed in 'equations.cc'. The attempted access was for variable with index " +
           std::to_string(var_index) + " .\n"));
  Assert(subset_attributes.at(var_index).var_type == field_type,
         dealii::ExcMessage(
           "PRISMS-PF Error: Attempted access of a " + field_string +
           " variable at the index of a non-" + field_string +
           " variable. The attempted access was for variable with index " +
           std::to_string(var_index) + " .\n"));
#endif
}

// The methods to set the residual terms
template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_scalar_value_term_RHS(
  unsigned int global_variable_index,
  T            val)
{
  scalar_vars_map[global_variable_index]->submit_value(val, q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_scalar_gradient_term_RHS(
  unsigned int              global_variable_index,
  dealii::Tensor<1, dim, T> grad)
{
  scalar_vars_map[global_variable_index]->submit_gradient(grad, q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_vector_value_term_RHS(
  unsigned int              global_variable_index,
  dealii::Tensor<1, dim, T> val)
{
  vector_vars_map[global_variable_index]->submit_value(val, q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_vector_gradient_term_RHS(
  unsigned int              global_variable_index,
  dealii::Tensor<2, dim, T> grad)
{
  vector_vars_map[global_variable_index]->submit_gradient(grad, q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_scalar_value_term_LHS(
  unsigned int global_variable_index,
  T            val)
{
  scalar_change_in_vars_map[global_variable_index]->submit_value(val, q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_scalar_gradient_term_LHS(
  unsigned int              global_variable_index,
  dealii::Tensor<1, dim, T> grad)
{
  scalar_change_in_vars_map[global_variable_index]->submit_gradient(grad, q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_vector_value_term_LHS(
  unsigned int              global_variable_index,
  dealii::Tensor<1, dim, T> val)
{
  vector_change_in_vars_map[global_variable_index]->submit_value(val, q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_vector_gradient_term_LHS(
  unsigned int              global_variable_index,
  dealii::Tensor<2, dim, T> grad)
{
  vector_change_in_vars_map[global_variable_index]->submit_gradient(grad, q_point);
}

template class variableContainer<2, 1, dealii::VectorizedArray<double>>;
template class variableContainer<3, 1, dealii::VectorizedArray<double>>;

template class variableContainer<2, 2, dealii::VectorizedArray<double>>;
template class variableContainer<3, 2, dealii::VectorizedArray<double>>;

template class variableContainer<2, 3, dealii::VectorizedArray<double>>;
template class variableContainer<3, 3, dealii::VectorizedArray<double>>;

template class variableContainer<2, 4, dealii::VectorizedArray<double>>;
template class variableContainer<3, 4, dealii::VectorizedArray<double>>;

template class variableContainer<2, 5, dealii::VectorizedArray<double>>;
template class variableContainer<3, 5, dealii::VectorizedArray<double>>;

template class variableContainer<2, 6, dealii::VectorizedArray<double>>;
template class variableContainer<3, 6, dealii::VectorizedArray<double>>;