#include <deal.II/matrix_free/evaluation_flags.h>

#include <config.h>
#include <core/exceptions.h>
#include <core/type_enums.h>
#include <core/variable_attributes.h>
#include <core/variable_container.h>
#include <string>

template <int dim, int degree, typename number>
variableContainer<dim, degree, number>::variableContainer(
  const dealii::MatrixFree<dim, number> &data,
  const AttributesList                  &_subset_attributes,
  const solveType                       &_solve_type)
  : subset_attributes(_subset_attributes)
  , solve_type(_solve_type)
{
  auto construct_map =
    [&](const std::unordered_map<std::pair<uint, dependencyType>, EvalFlags, pairHash>
          &eval_flag_set)
  {
    for (const auto &[pair, flags] : eval_flag_set)
      {
        const uint &dependency_index = pair.first;
        if (subset_attributes.at(dependency_index).field_type == fieldType::SCALAR)
          {
            scalar_vars_map[dependency_index].emplace(
              pair.second,
              std::make_unique<scalar_FEEval>(data, dependency_index));
          }
        else
          {
            vector_vars_map[dependency_index].emplace(
              pair.second,
              std::make_unique<vector_FEEval>(data, dependency_index));
          }
      }
  };

  // Loop through the variable attributes
  for (const auto &[index, variable] : subset_attributes)
    {
      if (solve_type == solveType::NONEXPLICIT_LHS)
        {
          construct_map(variable.eval_flag_set_LHS);
        }
      else
        {
          construct_map(variable.eval_flag_set_RHS);
        }
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(variableContainer &,
                           const dealii::Point<dim, dealii::VectorizedArray<number>> &)>
                                                     &func,
  dealii::LinearAlgebra::distributed::Vector<number> &dst,
  const LinearAlgebra::distributed::Vector<number>   &src,
  const std::pair<uint, uint>                        &cell_range)
{
  for (uint cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);

      for (uint q = 0; q < get_n_q_points(); ++q)
        {
          // Set the quadrature point
          q_point = q;

          // Grab the quadrature point location
          Point<dim, VectorizedArray<number>> q_point_loc = get_q_point_location();

          // Calculate the residuals
          func(*this, q_point_loc);
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval_local_diagonal(
  const std::function<void(variableContainer &,
                           const dealii::Point<dim, dealii::VectorizedArray<number>> &)>
                                                     &func,
  dealii::LinearAlgebra::distributed::Vector<number> &dst,
  const std::pair<uint, uint>                        &cell_range,
  const uint                                         &global_var_index)
{
  const auto &variable = subset_attributes.at(global_var_index);

  AssertThrow(variable.field_type != fieldType::VECTOR,
              FeatureNotImplemented("Vector multigrid"));

  auto *scalar_FEEval_ptr =
    scalar_vars_map.at(global_var_index).at(dependencyType::CHANGE).get();

  n_dofs_per_cell = scalar_FEEval_ptr->dofs_per_cell;
  diagonal = std::make_unique<dealii::AlignedVector<dealii::VectorizedArray<number>>>(
    n_dofs_per_cell);

  for (uint cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      reinit(cell, global_var_index);

      for (uint i = 0; i < n_dofs_per_cell; ++i)
        {
          for (uint j = 0; j < n_dofs_per_cell; ++j)
            {
              scalar_FEEval_ptr->submit_dof_value(dealii::VectorizedArray<number>(), j);
            }
          scalar_FEEval_ptr->submit_dof_value(dealii::make_vectorized_array<number>(1.0),
                                              i);

          eval(global_var_index);

          for (uint q = 0; q < get_n_q_points(); ++q)
            {
              // Set the quadrature point
              q_point = q;

              // Grab the quadrature point location
              dealii::Point<dim, dealii::VectorizedArray<number>> q_point_loc =
                get_q_point_location();

              // Calculate the residuals
              func(*this, q_point_loc);
            }

          integrate(global_var_index);
          (*diagonal)[i] = scalar_FEEval_ptr->get_dof_value(i);
        }

      for (uint i = 0; i < n_dofs_per_cell; ++i)
        {
          scalar_FEEval_ptr->submit_dof_value((*diagonal)[i], i);
        }
      scalar_FEEval_ptr->distribute_local_to_global(dst);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::access_valid(
  [[maybe_unused]] const uint           &global_variable_index,
  [[maybe_unused]] const dependencyType &dependency_type,
  [[maybe_unused]] const EvalFlags      &flag) const
{
  for ([[maybe_unused]] const auto &[index, variable] : subset_attributes)
    {
      if (solve_type == solveType::NONEXPLICIT_LHS)
        {
          Assert(
            variable.eval_flag_set_LHS.find(
              std::pair<uint, dependencyType>(global_variable_index, dependency_type)) !=
              variable.eval_flag_set_LHS.end(),
            dealii::ExcMessage("PRISMS-PF Error: Attempted access of a variable that was "
                               "not marked as needed in 'equations.cc'. The attempted "
                               "access was for variable with index " +
                               std::to_string(global_variable_index) +
                               "with the following dependency type" +
                               to_string(dependency_type) + "."));
        }
      else
        {
          Assert(
            variable.eval_flag_set_RHS.find(
              std::pair<uint, dependencyType>(global_variable_index, dependency_type)) !=
              variable.eval_flag_set_RHS.end(),
            dealii::ExcMessage("PRISMS-PF Error: Attempted access of a variable that was "
                               "not marked as needed in 'equations.cc'. The attempted "
                               "access was for variable with index " +
                               std::to_string(global_variable_index) +
                               "with the following dependency type" +
                               to_string(dependency_type) + "."));
        }
    }
}

template <int dim, int degree, typename number>
uint
variableContainer<dim, degree, number>::get_n_q_points() const
{
  if (!scalar_vars_map.empty())
    {
      const auto &first_scalar_FEEval = scalar_vars_map.begin()->second;

      Assert(!first_scalar_FEEval.empty(),
             dealii::ExcMessage(
               "PRISMS-PF Error: When trying to access the number of quadrature "
               "points, all FEEvaluation object containers were empty."));

      return first_scalar_FEEval.begin()->second->n_q_points;
    }
  else if (!vector_vars_map.empty())
    {
      const auto &first_vector_FEEval = vector_vars_map.begin()->second;

      Assert(!first_vector_FEEval.empty(),
             dealii::ExcMessage(
               "PRISMS-PF Error: When trying to access the number of quadrature "
               "points, all FEEvaluation object containers were empty."));

      return first_vector_FEEval.begin()->second->n_q_points;
    }

  Assert(false,
         dealii::ExcMessage(
           "PRISMS-PF Error: When trying to access the number of quadrature "
           "points, all FEEvaluation object containers were empty."));
  return 0;
}

template <int dim, int degree, typename number>
dealii::Point<dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_q_point_location() const
{
  if (!scalar_vars_map.empty())
    {
      const auto &first_scalar_FEEval = scalar_vars_map.begin()->second;

      Assert(!first_scalar_FEEval.empty(),
             dealii::ExcMessage(
               "PRISMS-PF Error: When trying to access the quadrature point "
               "location, all FEEvaluation object containers were empty."));

      return first_scalar_FEEval.begin()->second->quadrature_point(q_point);
    }
  else if (!vector_vars_map.empty())
    {
      const auto &first_vector_FEEval = vector_vars_map.begin()->second;

      Assert(!first_vector_FEEval.empty(),
             dealii::ExcMessage(
               "PRISMS-PF Error: When trying to access the quadrature point "
               "location, all FEEvaluation object containers were empty."));

      return first_vector_FEEval.begin()->second->quadrature_point(q_point);
    }

  Assert(false,
         dealii::ExcMessage("PRISMS-PF Error: When trying to access the quadrature point "
                            "location, all FEEvaluation object containers were empty."));

  return dealii::Point<dim, dealii::VectorizedArray<number>>();
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::reinit_and_eval(
  const dealii::LinearAlgebra::distributed::Vector<number> &src,
  uint                                                      cell)
{
  auto reinit_and_eval_map =
    [&](const std::unordered_map<std::pair<uint, dependencyType>, EvalFlags, pairHash>
          &eval_flag_set)
  {
    for (const auto &[pair, flags] : eval_flag_set)
      {
        const uint           &dependency_index = pair.first;
        const dependencyType &dependency_type  = pair.second;

        if (subset_attributes.at(dependency_index).field_type == fieldType::SCALAR)
          {
            auto *scalar_FEEval_ptr =
              scalar_vars_map.at(dependency_index).at(dependency_type).get();
            scalar_FEEval_ptr->reinit(cell);
            scalar_FEEval_ptr->read_dof_values_plain(src);
            scalar_FEEval_ptr->evaluate(flags);
          }
        else
          {
            auto *vector_FEEval_ptr =
              vector_vars_map.at(dependency_index).at(dependency_type).get();
            vector_FEEval_ptr->reinit(cell);
            vector_FEEval_ptr->read_dof_values_plain(src);
            vector_FEEval_ptr->evaluate(flags);
          }
      }
  };

  for (const auto &[var_index, variable] : subset_attributes)
    {
      if (solve_type == solveType::NONEXPLICIT_LHS)
        {
          reinit_and_eval_map(variable.eval_flag_set_LHS);
        }
      else
        {
          reinit_and_eval_map(variable.eval_flag_set_RHS);
        }
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::reinit(uint        cell,
                                               const uint &global_variable_index)
{
  auto reinit_map =
    [&](const std::unordered_map<std::pair<uint, dependencyType>, EvalFlags, pairHash>
          &eval_flag_set)
  {
    for (const auto &[pair, flags] : eval_flag_set)
      {
        const uint           &dependency_index = pair.first;
        const dependencyType &dependency_type  = pair.second;

        if (subset_attributes.at(dependency_index).field_type == fieldType::SCALAR)
          {
            auto *scalar_FEEval_ptr =
              scalar_vars_map.at(dependency_index).at(dependency_type).get();
            scalar_FEEval_ptr->reinit(cell);
          }
        else
          {
            auto *vector_FEEval_ptr =
              vector_vars_map.at(dependency_index).at(dependency_type).get();
            vector_FEEval_ptr->reinit(cell);
          }
      }
  };

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      reinit_map(subset_attributes.at(global_variable_index).eval_flag_set_LHS);
    }
  else
    {
      reinit_map(subset_attributes.at(global_variable_index).eval_flag_set_RHS);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval(const uint &global_variable_index)
{
  auto eval_map =
    [&](const std::unordered_map<std::pair<uint, dependencyType>, EvalFlags, pairHash>
          &eval_flag_set)
  {
    for (const auto &[pair, flags] : eval_flag_set)
      {
        const uint           &dependency_index = pair.first;
        const dependencyType &dependency_type  = pair.second;

        if (subset_attributes.at(dependency_index).field_type == fieldType::SCALAR)
          {
            auto *scalar_FEEval_ptr =
              scalar_vars_map.at(dependency_index).at(dependency_type).get();
            scalar_FEEval_ptr->evaluate(flags);
          }
        else
          {
            auto *vector_FEEval_ptr =
              vector_vars_map.at(dependency_index).at(dependency_type).get();
            vector_FEEval_ptr->evaluate(flags);
          }
      }
  };

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      eval_map(subset_attributes.at(global_variable_index).eval_flag_set_LHS);
    }
  else
    {
      eval_map(subset_attributes.at(global_variable_index).eval_flag_set_RHS);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::integrate(const uint &global_variable_index)
{
  const auto &variable = subset_attributes.at(global_variable_index);

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      if (variable.field_type == fieldType::SCALAR)
        {
          auto *scalar_FEEval_ptr =
            scalar_vars_map.at(global_variable_index).at(dependencyType::CHANGE).get();
          scalar_FEEval_ptr->integrate(variable.eval_flags_residual_LHS);
        }
      else
        {
          auto *vector_FEEval_ptr =
            vector_vars_map.at(global_variable_index).at(dependencyType::CHANGE).get();
          vector_FEEval_ptr->integrate(variable.eval_flags_residual_LHS);
        }
    }
  else
    {
      if (subset_attributes.at(global_variable_index).field_type == fieldType::SCALAR)
        {
          auto *scalar_FEEval_ptr =
            scalar_vars_map.at(global_variable_index).at(dependencyType::CHANGE).get();
          scalar_FEEval_ptr->integrate(variable.eval_flags_residual_RHS);
        }
      else
        {
          auto *vector_FEEval_ptr =
            vector_vars_map.at(global_variable_index).at(dependencyType::CHANGE).get();
          vector_FEEval_ptr->integrate(variable.eval_flags_residual_RHS);
        }
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::integrate_and_distribute(
  dealii::LinearAlgebra::distributed::Vector<number> &dst)
{
  auto integrate_and_distribute_map = [&](const EvalFlags      &residual_flag_set,
                                          const dependencyType &dependency_type,
                                          const uint           &residual_index)
  {
    if (subset_attributes.at(residual_index).field_type == fieldType::SCALAR)
      {
        auto *scalar_FEEval_ptr =
          scalar_vars_map.at(residual_index).at(dependency_type).get();
        scalar_FEEval_ptr->integrate_scatter(residual_flag_set, dst);
      }
    else
      {
        auto *vector_FEEval_ptr =
          vector_vars_map.at(residual_index).at(dependency_type).get();
        vector_FEEval_ptr->integrate_scatter(residual_flag_set, dst);
      }
  };

  for (const auto &[index, variable] : subset_attributes)
    {
      if (solve_type == solveType::NONEXPLICIT_LHS)
        {
          integrate_and_distribute_map(variable.eval_flags_residual_LHS,
                                       dependencyType::CHANGE,
                                       index);
        }
      else
        {
          integrate_and_distribute_map(variable.eval_flags_residual_RHS,
                                       dependencyType::NORMAL,
                                       index);
        }
    }
}

template <int dim, int degree, typename number>
dealii::VectorizedArray<number>
variableContainer<dim, degree, number>::get_scalar_value(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::values);
#endif

  return scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_value(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_scalar_gradient(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::gradients);
#endif

  return scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_gradient(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_scalar_hessian(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::hessians);
#endif

  return scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_hessian(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_scalar_hessian_diagonal(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::hessians);
#endif

  return scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_hessian_diagonal(q_point);
}

template <int dim, int degree, typename number>
dealii::VectorizedArray<number>
variableContainer<dim, degree, number>::get_scalar_laplacian(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::hessians);
#endif

  return scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_laplacian(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_vector_value(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::values);
#endif

  const auto &value =
    vector_vars_map.at(global_variable_index).at(dependency_type)->get_value(q_point);

  if constexpr (dim == 1)
    {
      // Wrap the value for consistency
      dealii::Tensor<1, dim, dealii::VectorizedArray<number>> wrapper;
      wrapper[0] = value;
      return wrapper;
    }
  else
    {
      // Return the value directly for dim > 1
      return value;
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_vector_gradient(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::gradients);
#endif

  const auto &grad =
    vector_vars_map.at(global_variable_index).at(dependency_type)->get_gradient(q_point);

  if constexpr (dim == 1)
    {
      // Wrap the value for consistency
      dealii::Tensor<2, dim, dealii::VectorizedArray<number>> wrapper;
      wrapper[0] = grad;
      return wrapper;
    }
  else
    {
      // Return the value directly for dim > 1
      return grad;
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<3, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_vector_hessian(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::hessians);
#endif

  const auto &hess =
    vector_vars_map.at(global_variable_index).at(dependency_type)->get_hessian(q_point);

  if constexpr (dim == 1)
    {
      // Wrap the value for consistency
      dealii::Tensor<3, dim, dealii::VectorizedArray<number>> wrapper;
      wrapper[0] = hess;
      return wrapper;
    }
  else
    {
      // Return the value directly for dim > 1
      return hess;
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_vector_hessian_diagonal(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::hessians);
#endif

  const auto &hess_diag = vector_vars_map.at(global_variable_index)
                            .at(dependency_type)
                            ->get_hessian_diagonal(q_point);

  if constexpr (dim == 1)
    {
      // Wrap the value for consistency
      dealii::Tensor<2, dim, dealii::VectorizedArray<number>> wrapper;
      wrapper[0] = hess_diag;
      return wrapper;
    }
  else
    {
      // Return the value directly for dim > 1
      return hess_diag;
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_vector_laplacian(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::hessians);
#endif

  const auto &lap =
    vector_vars_map.at(global_variable_index).at(dependency_type)->get_laplacian(q_point);

  if constexpr (dim == 1)
    {
      // Wrap the value for consistency
      dealii::Tensor<1, dim, dealii::VectorizedArray<number>> wrapper;
      wrapper[0] = lap;
      return wrapper;
    }
  else
    {
      // Return the value directly for dim > 1
      return lap;
    }
}

template <int dim, int degree, typename number>
dealii::VectorizedArray<number>
variableContainer<dim, degree, number>::get_vector_divergence(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::gradients);
#endif

  return vector_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_divergence(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_vector_symmetric_gradient(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::gradients);
#endif

  return vector_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_symmetric_gradient(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1, (dim == 2 ? 1 : dim), dealii::VectorizedArray<number>>
variableContainer<dim, degree, number>::get_vector_curl(
  uint           global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index, dependency_type, EvalFlags::gradients);
#endif

  if constexpr (dim == 1)
    {
      Assert(false, dealii::ExcMessage("PRISMS-PF Error: Curl is nonsensical for 1D."));
      return dealii::Tensor<1, (dim == 2 ? 1 : dim), dealii::VectorizedArray<number>> {};
    }
  else
    {
      // Return the value directly for dim > 1
      return vector_vars_map.at(global_variable_index)
        .at(dependency_type)
        ->get_curl(q_point);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_scalar_value_term(
  uint                            global_variable_index,
  dealii::VectorizedArray<number> val,
  dependencyType                  dependency_type)
{
#ifdef DEBUG
  Assert(
    dependency_type == NORMAL || solve_type == solveType::NONEXPLICIT_LHS,
    dealii::ExcMessage(
      "PRISMS-PF Error: RHS residuals are only allowed to submit normal value terms."));
  Assert(
    dependency_type == CHANGE || solve_type != solveType::NONEXPLICIT_LHS,
    dealii::ExcMessage(
      "PRISMS-PF Error: LHS residuals are only allowed to submit change value terms."));
#endif

  scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->submit_value(val, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_scalar_gradient_term(
  uint                                                    global_variable_index,
  dealii::Tensor<1, dim, dealii::VectorizedArray<number>> grad,
  dependencyType                                          dependency_type)
{
#ifdef DEBUG
  Assert(dependency_type == NORMAL || solve_type == solveType::NONEXPLICIT_LHS,
         dealii::ExcMessage("PRISMS-PF Error: RHS residuals are only allowed to submit "
                            "normal gradient terms."));
  Assert(dependency_type == CHANGE || solve_type != solveType::NONEXPLICIT_LHS,
         dealii::ExcMessage("PRISMS-PF Error: LHS residuals are only allowed to submit "
                            "change gradient terms."));
#endif

  scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->submit_gradient(grad, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_vector_value_term(
  uint                                                    global_variable_index,
  dealii::Tensor<1, dim, dealii::VectorizedArray<number>> val,
  dependencyType                                          dependency_type)
{
#ifdef DEBUG
  Assert(
    dependency_type == NORMAL || solve_type == solveType::NONEXPLICIT_LHS,
    dealii::ExcMessage(
      "PRISMS-PF Error: RHS residuals are only allowed to submit normal value terms."));
  Assert(
    dependency_type == CHANGE || solve_type != solveType::NONEXPLICIT_LHS,
    dealii::ExcMessage(
      "PRISMS-PF Error: LHS residuals are only allowed to submit change value terms."));
#endif

  vector_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->submit_value(val, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_vector_gradient_term(
  uint                                                    global_variable_index,
  dealii::Tensor<2, dim, dealii::VectorizedArray<number>> grad,
  dependencyType                                          dependency_type)
{
#ifdef DEBUG
  Assert(dependency_type == NORMAL || solve_type == solveType::NONEXPLICIT_LHS,
         dealii::ExcMessage("PRISMS-PF Error: RHS residuals are only allowed to submit "
                            "normal gradient terms."));
  Assert(dependency_type == CHANGE || solve_type != solveType::NONEXPLICIT_LHS,
         dealii::ExcMessage("PRISMS-PF Error: LHS residuals are only allowed to submit "
                            "change gradient terms."));
#endif

  vector_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->submit_gradient(grad, q_point);
}

INSTANTIATE_TRI_TEMPLATE(variableContainer)