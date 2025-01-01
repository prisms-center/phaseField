#include <deal.II/base/exceptions.h>

#include <core/variableContainer.h>
#include <string>

template <int dim, int degree, typename number>
variableContainer<dim, degree, number>::variableContainer(
  const dealii::MatrixFree<dim, number> &data,
  const std::vector<variable_info>      &_variable_info_list,
  const std::vector<variable_info>      &_old_variable_info_list,
  const std::vector<variable_info>      &_change_variable_info_list)
  : variable_info_list(_variable_info_list)
  , old_variable_info_list(_old_variable_info_list)
  , change_variable_info_list(_change_variable_info_list)
  , num_var(variable_info_list.size())
{
  for (unsigned int i = 0; i < num_var; i++)
    {
      const auto &var_info        = variable_info_list[i];
      const auto &var_change_info = change_variable_info_list[i];

      if (var_info.var_needed)
        {
          const unsigned int var_index = var_info.global_var_index;

          if (var_info.is_scalar)
            {
              scalar_vars_map.emplace(var_index,
                                      std::make_unique<scalar_FEEval>(data, i));
            }
          else
            {
              vector_vars_map.emplace(var_index,
                                      std::make_unique<vector_FEEval>(data, i));
            }
        }

      if (var_change_info.var_needed)
        {
          const unsigned int var_index = var_change_info.global_var_index;

          if (var_change_info.is_scalar)
            {
              scalar_change_in_vars_map.emplace(var_index,
                                                std::make_unique<scalar_FEEval>(data, i));
            }
          else
            {
              vector_change_in_vars_map.emplace(var_index,
                                                std::make_unique<vector_FEEval>(data, i));
            }
        }
    }
}

template <int dim, int degree, typename number>
variableContainer<dim, degree, number>::variableContainer(
  const dealii::MatrixFree<dim, number> &data,
  const std::vector<variable_info>      &_variable_info_list,
  const std::vector<variable_info>      &_old_variable_info_list)
  : variable_info_list(_variable_info_list)
  , old_variable_info_list(_old_variable_info_list)
  , num_var(variable_info_list.size())
{
  for (unsigned int i = 0; i < num_var; i++)
    {
      const auto &var_info = variable_info_list[i];

      if (!var_info.var_needed)
        {
          continue;
        }

      const unsigned int var_index = var_info.global_var_index;

      if (var_info.is_scalar)
        {
          scalar_vars_map.emplace(var_index, std::make_unique<scalar_FEEval>(data, i));
        }
      else
        {
          vector_vars_map.emplace(var_index, std::make_unique<vector_FEEval>(data, i));
        }
    }
}

// Variant of the constructor where it reads from a fixed index of "data", used
// for post-processing
template <int dim, int degree, typename number>
variableContainer<dim, degree, number>::variableContainer(
  const dealii::MatrixFree<dim, number> &data,
  const std::vector<variable_info>      &_variable_info_list,
  const unsigned int                    &fixed_index)
  : variable_info_list(_variable_info_list)
  , num_var(variable_info_list.size())
{
  for (unsigned int i = 0; i < num_var; i++)
    {
      const auto &var_info = variable_info_list[i];

      if (!var_info.var_needed)
        {
          continue;
        }

      const unsigned int var_index = var_info.global_var_index;

      if (var_info.is_scalar)
        {
          scalar_vars_map.emplace(var_index,
                                  std::make_unique<scalar_FEEval>(data, fixed_index));
        }
      else
        {
          vector_vars_map.emplace(var_index,
                                  std::make_unique<vector_FEEval>(data, fixed_index));
        }
    }
}

template <int dim, int degree, typename number>
unsigned int
variableContainer<dim, degree, number>::get_num_q_points() const
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

template <int dim, int degree, typename number>
dealii::Point<dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_q_point_location() const
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

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::reinit_and_eval(
  const std::vector<dealii::LinearAlgebra::distributed::Vector<number> *> &src,
  unsigned int                                                             cell)
{
  for (unsigned int i = 0; i < num_var; i++)
    {
      const auto &var_info = variable_info_list[i];

      if (!var_info.var_needed)
        {
          continue;
        }

      const unsigned int var_index = var_info.global_var_index;

      if (var_info.is_scalar)
        {
          auto *scalar_FEEval_ptr = scalar_vars_map[var_index].get();
          scalar_FEEval_ptr->reinit(cell);
          scalar_FEEval_ptr->read_dof_values(*src[i]);
          scalar_FEEval_ptr->evaluate(var_info.evaluation_flags);
        }
      else
        {
          auto *vector_FEEval_ptr = vector_vars_map[var_index].get();
          vector_FEEval_ptr->reinit(cell);
          vector_FEEval_ptr->read_dof_values(*src[i]);
          vector_FEEval_ptr->evaluate(var_info.evaluation_flags);
        }
    }
}

/**
 * This is specialized for the LHS where a change in solution is needed. The RHS method
 * takes the src as a vector of LinearAlgebra::distributed::Vector<double>.
 */
template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::reinit_and_eval_change_in_solution(
  const dealii::LinearAlgebra::distributed::Vector<number> &src,
  unsigned int                                              cell,
  unsigned int                                              var_being_solved)
{
  if (change_variable_info_list[var_being_solved].is_scalar)
    {
      auto *scalar_FEEval_ptr = scalar_change_in_vars_map[var_being_solved].get();
      scalar_FEEval_ptr->reinit(cell);
      scalar_FEEval_ptr->read_dof_values(src);
      scalar_FEEval_ptr->evaluate(
        change_variable_info_list[var_being_solved].evaluation_flags);
    }
  else
    {
      auto *vector_FEEval_ptr = vector_change_in_vars_map[var_being_solved].get();
      vector_FEEval_ptr->reinit(cell);
      vector_FEEval_ptr->read_dof_values(src);
      vector_FEEval_ptr->evaluate(
        change_variable_info_list[var_being_solved].evaluation_flags);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::reinit(unsigned int cell)
{
  for (unsigned int i = 0; i < num_var; i++)
    {
      const auto &var_info = variable_info_list[i];

      if (!var_info.var_needed)
        {
          continue;
        }

      const unsigned int var_index = var_info.global_var_index;

      if (var_info.is_scalar)
        {
          scalar_vars_map[var_index]->reinit(cell);
        }
      else
        {
          vector_vars_map[var_index]->reinit(cell);
        }
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::integrate_and_distribute(
  std::vector<dealii::LinearAlgebra::distributed::Vector<number> *> &dst)
{
  for (unsigned int i = 0; i < num_var; i++)
    {
      const auto &var_info = variable_info_list[i];

      if (var_info.residual_flags & dealii::EvaluationFlags::nothing)
        {
          continue;
        }

      const unsigned int var_index = var_info.global_var_index;

      if (var_info.is_scalar)
        {
          auto *scalar_FEEval_ptr = scalar_vars_map[var_index].get();
          scalar_FEEval_ptr->integrate_scatter(var_info.residual_flags, *dst[i]);
        }
      else
        {
          auto *vector_FEEval_ptr = vector_vars_map[var_index].get();
          vector_FEEval_ptr->integrate_scatter(var_info.residual_flags, *dst[i]);
        }
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::integrate_and_distribute_change_in_solution_LHS(
  dealii::LinearAlgebra::distributed::Vector<number> &dst,
  const unsigned int                                  var_being_solved)
{
  if (change_variable_info_list[var_being_solved].is_scalar)
    {
      auto *scalar_FEEval_ptr = scalar_change_in_vars_map[var_being_solved].get();
      scalar_FEEval_ptr->integrate_scatter(
        change_variable_info_list[var_being_solved].residual_flags,
        dst);
    }
  else
    {
      auto *vector_FEEval_ptr = vector_change_in_vars_map[var_being_solved].get();
      vector_FEEval_ptr->integrate_scatter(
        change_variable_info_list[var_being_solved].residual_flags,
        dst);
    }
}

// Need to add index checking to these functions so that an error is thrown if
// the index wasn't set
template <int dim, int degree, typename number>
dealii::VectorizedArray<number>
variableContainer<dim, degree, number>::get_scalar_value(
  unsigned int global_variable_index) const
{
  Assert(!(variable_info_list[global_variable_index].evaluation_flags &
           dealii::EvaluationFlags::values),
         dealii::ExcMessage(
           "PRISMS-PF Error: Attempted access of a variable value that was not marked as "
           "needed in 'equations.cc'. The attempted access was for variable with index " +
           std::to_string(global_variable_index) + "."));

  return scalar_vars_map.at(global_variable_index)->get_value(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_scalar_gradient(
  unsigned int global_variable_index) const
{
  Assert(
    !(variable_info_list[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::gradients),
    dealii::ExcMessage(
      "PRISMS-PF Error: Attempted access of a variable gradient that was not marked as "
      "needed in 'equations.cc'. The attempted access was for variable with index " +
      std::to_string(global_variable_index) + "."));

  return scalar_vars_map.at(global_variable_index)->get_gradient(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_scalar_hessian(
  unsigned int global_variable_index) const
{
  Assert(
    !(variable_info_list[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::hessians),
    dealii::ExcMessage(
      "PRISMS-PF Error: Attempted access of a variable hessian that was not marked as "
      "needed in 'equations.cc'. The attempted access was for variable with index " +
      std::to_string(global_variable_index) + "."));

  return scalar_vars_map.at(global_variable_index)->get_hessian(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_vector_value(
  unsigned int global_variable_index) const
{
  Assert(!(variable_info_list[global_variable_index].evaluation_flags &
           dealii::EvaluationFlags::values),
         dealii::ExcMessage(
           "PRISMS-PF Error: Attempted access of a variable value that was not marked as "
           "needed in 'equations.cc'. The attempted access was for variable with index " +
           std::to_string(global_variable_index) + "."));

  return vector_vars_map.at(global_variable_index)->get_value(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_vector_gradient(
  unsigned int global_variable_index) const
{
  Assert(
    !(variable_info_list[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::gradients),
    dealii::ExcMessage(
      "PRISMS-PF Error: Attempted access of a variable gradient that was not marked as "
      "needed in 'equations.cc'. The attempted access was for variable with index " +
      std::to_string(global_variable_index) + "."));

  return vector_vars_map.at(global_variable_index)->get_gradient(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<3, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_vector_hessian(
  unsigned int global_variable_index) const
{
  Assert(
    !(variable_info_list[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::hessians),
    dealii::ExcMessage(
      "PRISMS-PF Error: Attempted access of a variable hessian that was not marked as "
      "needed in 'equations.cc'. The attempted access was for variable with index " +
      std::to_string(global_variable_index) + "."));

  return vector_vars_map.at(global_variable_index)->get_hessian(q_point);
}

// Need to add index checking to these functions so that an error is thrown if
// the index wasn't set
template <int dim, int degree, typename number>
dealii::VectorizedArray<number>
variableContainer<dim, degree, number>::get_change_in_scalar_value(
  unsigned int global_variable_index) const
{
  Assert(!(change_variable_info_list[global_variable_index].evaluation_flags &
           dealii::EvaluationFlags::values),
         dealii::ExcMessage(
           "PRISMS-PF Error: Attempted access of a change in variable value that was not "
           "marked as needed in 'equations.cc'. The attempted access was for variable "
           "with index " +
           std::to_string(global_variable_index) + "."));

  return scalar_change_in_vars_map.at(global_variable_index)->get_value(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_change_in_scalar_gradient(
  unsigned int global_variable_index) const
{
  Assert(
    !(change_variable_info_list[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::gradients),
    dealii::ExcMessage(
      "PRISMS-PF Error: Attempted access of a change in variable gradient that was not "
      "marked as needed in 'equations.cc'. The attempted access was for variable "
      "with index " +
      std::to_string(global_variable_index) + "."));

  return scalar_change_in_vars_map.at(global_variable_index)->get_gradient(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_change_in_scalar_hessian(
  unsigned int global_variable_index) const
{
  Assert(
    !(change_variable_info_list[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::hessians),
    dealii::ExcMessage(
      "PRISMS-PF Error: Attempted access of a change in variable hessian that was not "
      "marked as needed in 'equations.cc'. The attempted access was for variable "
      "with index " +
      std::to_string(global_variable_index) + "."));

  return scalar_change_in_vars_map.at(global_variable_index)->get_hessian(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_change_in_vector_value(
  unsigned int global_variable_index) const
{
  Assert(!(change_variable_info_list[global_variable_index].evaluation_flags &
           dealii::EvaluationFlags::values),
         dealii::ExcMessage(
           "PRISMS-PF Error: Attempted access of a change in variable value that was not "
           "marked as needed in 'equations.cc'. The attempted access was for variable "
           "with index " +
           std::to_string(global_variable_index) + "."));

  return vector_change_in_vars_map.at(global_variable_index)->get_value(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_change_in_vector_gradient(
  unsigned int global_variable_index) const
{
  Assert(
    !(change_variable_info_list[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::gradients),
    dealii::ExcMessage(
      "PRISMS-PF Error: Attempted access of a change in variable gradient that was not "
      "marked as needed in 'equations.cc'. The attempted access was for variable "
      "with index " +
      std::to_string(global_variable_index) + "."));

  return vector_change_in_vars_map.at(global_variable_index)->get_gradient(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<3, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_change_in_vector_hessian(
  unsigned int global_variable_index) const
{
  Assert(
    !(change_variable_info_list[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::hessians),
    dealii::ExcMessage(
      "PRISMS-PF Error: Attempted access of a change in variable hessian that was not "
      "marked as needed in 'equations.cc'. The attempted access was for variable "
      "with index " +
      std::to_string(global_variable_index) + "."));

  return vector_change_in_vars_map.at(global_variable_index)->get_hessian(q_point);
}

// The methods to set the residual terms
template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_scalar_value_term_RHS(
  unsigned int                    global_variable_index,
  dealii::VectorizedArray<number> val)
{
  scalar_vars_map[global_variable_index]->submit_value(val, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_scalar_gradient_term_RHS(
  unsigned int                                            global_variable_index,
  dealii::Tensor<1, dim, dealii::VectorizedArray<number>> grad)
{
  scalar_vars_map[global_variable_index]->submit_gradient(grad, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_vector_value_term_RHS(
  unsigned int                                            global_variable_index,
  dealii::Tensor<1, dim, dealii::VectorizedArray<number>> val)
{
  vector_vars_map[global_variable_index]->submit_value(val, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_vector_gradient_term_RHS(
  unsigned int                                            global_variable_index,
  dealii::Tensor<2, dim, dealii::VectorizedArray<number>> grad)
{
  vector_vars_map[global_variable_index]->submit_gradient(grad, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_scalar_value_term_LHS(
  unsigned int                    global_variable_index,
  dealii::VectorizedArray<number> val)
{
  scalar_change_in_vars_map[global_variable_index]->submit_value(val, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_scalar_gradient_term_LHS(
  unsigned int                                            global_variable_index,
  dealii::Tensor<1, dim, dealii::VectorizedArray<number>> grad)
{
  scalar_change_in_vars_map[global_variable_index]->submit_gradient(grad, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_vector_value_term_LHS(
  unsigned int                                            global_variable_index,
  dealii::Tensor<1, dim, dealii::VectorizedArray<number>> val)
{
  vector_change_in_vars_map[global_variable_index]->submit_value(val, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_vector_gradient_term_LHS(
  unsigned int                                            global_variable_index,
  dealii::Tensor<2, dim, dealii::VectorizedArray<number>> grad)
{
  vector_change_in_vars_map[global_variable_index]->submit_gradient(grad, q_point);
}
