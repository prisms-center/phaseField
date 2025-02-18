#include <deal.II/base/point.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/config.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/core/variable_container.h>

#include <string>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, int degree, typename number>
variableContainer<dim, degree, number>::variableContainer(
  const dealii::MatrixFree<dim, number>            &data,
  const std::map<unsigned int, variableAttributes> &_subset_attributes,
  const std::unordered_map<std::pair<unsigned int, dependencyType>,
                           unsigned int,
                           pairHash>               &_global_to_local_solution,
  const solveType                                  &_solve_type)
  : subset_attributes(_subset_attributes)
  , global_to_local_solution(_global_to_local_solution)
  , solve_type(_solve_type)
{
  auto construct_map =
    [&](const std::map<unsigned int, std::map<dependencyType, fieldType>> &dependency_set)
  {
    for (const auto &[dependency_index, map] : dependency_set)
      {
        for (const auto &[dependency_type, field_type] : map)
          {
            if (field_type == fieldType::SCALAR)
              {
                scalar_vars_map[dependency_index].emplace(
                  dependency_type,
                  std::make_unique<scalar_FEEval>(data, dependency_index));
              }
            else
              {
                vector_vars_map[dependency_index].emplace(
                  dependency_type,
                  std::make_unique<vector_FEEval>(data, dependency_index));
              }
          }
      }
  };

  // For explicit solves we have already flattened the dependencies
  if (solve_type == solveType::EXPLICIT_RHS || solve_type == solveType::POSTPROCESS)
    {
      construct_map(subset_attributes.begin()->second.dependency_set_RHS);
      return;
    }

  // TODO: Add stuff for cononlinear solves

  // Loop through the variable attributes for nonexplicit solves
  Assert(subset_attributes.size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      construct_map(subset_attributes.begin()->second.dependency_set_LHS);
    }
  else
    {
      construct_map(subset_attributes.begin()->second.dependency_set_RHS);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
                                              &func,
  std::vector<VectorType *>                   &dst,
  const std::vector<VectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);

      for (unsigned int q = 0; q < get_n_q_points(); ++q)
        {
          // Set the quadrature point
          q_point = q;

          // Grab the quadrature point location
          dealii::Point<dim, size_type> q_point_loc = get_q_point_location();

          // Calculate the residuals
          func(*this, q_point_loc);
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
                                              &func,
  VectorType                                  &dst,
  const std::vector<VectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);

      for (unsigned int q = 0; q < get_n_q_points(); ++q)
        {
          // Set the quadrature point
          q_point = q;

          // Grab the quadrature point location
          dealii::Point<dim, size_type> q_point_loc = get_q_point_location();

          // Calculate the residuals
          func(*this, q_point_loc);
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
                                              &func,
  VectorType                                  &dst,
  const VectorType                            &src,
  const std::vector<VectorType *>             &src_subset,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);
      reinit_and_eval(src_subset, cell);

      for (unsigned int q = 0; q < get_n_q_points(); ++q)
        {
          // Set the quadrature point
          q_point = q;

          // Grab the quadrature point location
          dealii::Point<dim, size_type> q_point_loc = get_q_point_location();

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
  const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
                                              &func,
  VectorType                                  &dst,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  Assert(subset_attributes.size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  const auto &global_var_index = subset_attributes.begin()->first;

  AssertThrow(subset_attributes.begin()->second.field_type != fieldType::VECTOR,
              FeatureNotImplemented("Vector multigrid"));

  auto *scalar_FEEval_ptr =
    scalar_vars_map.at(global_var_index).at(dependencyType::CHANGE).get();

  n_dofs_per_cell = scalar_FEEval_ptr->dofs_per_cell;
  diagonal        = std::make_unique<dealii::AlignedVector<size_type>>(n_dofs_per_cell);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      reinit(cell, global_var_index);

      for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
            {
              scalar_FEEval_ptr->submit_dof_value(size_type(), j);
            }
          scalar_FEEval_ptr->submit_dof_value(dealii::make_vectorized_array<number>(1.0),
                                              i);

          eval(global_var_index);

          for (unsigned int q = 0; q < get_n_q_points(); ++q)
            {
              // Set the quadrature point
              q_point = q;

              // Grab the quadrature point location
              dealii::Point<dim, size_type> q_point_loc = get_q_point_location();

              // Calculate the residuals
              func(*this, q_point_loc);
            }

          integrate(global_var_index);
          (*diagonal)[i] = scalar_FEEval_ptr->get_dof_value(i);
        }

      for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
        {
          scalar_FEEval_ptr->submit_dof_value((*diagonal)[i], i);
        }
      scalar_FEEval_ptr->distribute_local_to_global(dst);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::scalar_FEEval_exists(
  [[maybe_unused]] const unsigned int   &dependency_index,
  [[maybe_unused]] const dependencyType &dependency_type) const
{
  Assert(scalar_vars_map.find(dependency_index) != scalar_vars_map.end(),
         dealii::ExcMessage(
           "The scalar FEEvaluation object does not exist for global index = " +
           std::to_string(dependency_index)));
  Assert(scalar_vars_map.at(dependency_index).find(dependency_type) !=
           scalar_vars_map.at(dependency_index).end(),
         dealii::ExcMessage("The scalar FEEvaluation object with global index = " +
                            std::to_string(dependency_index) +
                            " does not exist for type = " + to_string(dependency_type)));
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::vector_FEEval_exists(
  [[maybe_unused]] const unsigned int   &dependency_index,
  [[maybe_unused]] const dependencyType &dependency_type) const
{
  Assert(vector_vars_map.find(dependency_index) != vector_vars_map.end(),
         dealii::ExcMessage(
           "The vector FEEvaluation object does not exist for global index = " +
           std::to_string(dependency_index)));
  Assert(vector_vars_map.at(dependency_index).find(dependency_type) !=
           vector_vars_map.at(dependency_index).end(),
         dealii::ExcMessage("The vector FEEvaluation object with global index = " +
                            std::to_string(dependency_index) +
                            " does not exist for type = " + to_string(dependency_type)));
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::access_valid(
  [[maybe_unused]] const unsigned int                             &dependency_index,
  [[maybe_unused]] const dependencyType                           &dependency_type,
  [[maybe_unused]] const dealii::EvaluationFlags::EvaluationFlags &flag) const
{
  for ([[maybe_unused]] const auto &[index, variable] : subset_attributes)
    {
      if (solve_type == solveType::NONEXPLICIT_LHS)
        {
          Assert(variable.eval_flag_set_LHS.find(
                   std::make_pair(dependency_index, dependency_type)) !=
                   variable.eval_flag_set_LHS.end(),
                 DependencyNotFound(dependency_index, to_string(dependency_type)));
        }
      else
        {
          Assert(variable.eval_flag_set_RHS.find(
                   std::make_pair(dependency_index, dependency_type)) !=
                   variable.eval_flag_set_RHS.end(),
                 DependencyNotFound(dependency_index, to_string(dependency_type)));
        }
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::submission_valid(
  [[maybe_unused]] const dependencyType &dependency_type) const
{
  Assert(dependency_type == NORMAL || solve_type == solveType::NONEXPLICIT_LHS,
         dealii::ExcMessage("PRISMS-PF Error: RHS residuals are only allowed to submit "
                            "normal gradient terms."));
  Assert(dependency_type == CHANGE || solve_type != solveType::NONEXPLICIT_LHS,
         dealii::ExcMessage("PRISMS-PF Error: LHS residuals are only allowed to submit "
                            "change gradient terms."));
}

template <int dim, int degree, typename number>
unsigned int
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
dealii::Point<dim, typename variableContainer<dim, degree, number>::size_type>
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

  return dealii::Point<dim, size_type>();
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::reinit_and_eval(
  const std::vector<VectorType *> &src,
  unsigned int                     cell)
{
  // Reinit and eval values for the given dependency set. Note the dependency set includes
  // the variable we're evaluating, which may or may not be an actually dependency. For
  // this reason, I selectively read dofs and evaluate the flags.
  auto reinit_and_eval_map =
    [&](const std::unordered_map<std::pair<unsigned int, dependencyType>,
                                 dealii::EvaluationFlags::EvaluationFlags,
                                 pairHash>                                &eval_flag_set,
        const std::map<unsigned int, std::map<dependencyType, fieldType>> &dependency_set)
  {
    for (const auto &[dependency_index, map] : dependency_set)
      {
        for (const auto &[dependency_type, field_type] : map)
          {
            if (dependency_type == dependencyType::CHANGE)
              {
                continue;
              }

            const auto &pair = std::make_pair(dependency_index, dependency_type);

            Assert(global_to_local_solution.find(pair) != global_to_local_solution.end(),
                   dealii::ExcMessage(
                     "The global to local mapping does not exists for global index = " +
                     std::to_string(dependency_index) +
                     "  and type = " + to_string(dependency_type)));

            if (field_type == fieldType::SCALAR)
              {
                scalar_FEEval_exists(dependency_index, dependency_type);

                auto *scalar_FEEval_ptr =
                  scalar_vars_map.at(dependency_index).at(dependency_type).get();
                scalar_FEEval_ptr->reinit(cell);

                if (eval_flag_set.find(pair) != eval_flag_set.end())
                  {
                    const unsigned int &local_index = global_to_local_solution.at(pair);

                    Assert(src.size() > local_index,
                           dealii::ExcMessage(
                             "The provided src vector's size is below the given local "
                             "index = " +
                             std::to_string(local_index) +
                             " for global index = " + std::to_string(dependency_index) +
                             "  and type = " + to_string(dependency_type)));

                    scalar_FEEval_ptr->read_dof_values_plain(*(src.at(local_index)));
                    scalar_FEEval_ptr->evaluate(eval_flag_set.at(pair));
                  }
              }
            else
              {
                vector_FEEval_exists(dependency_index, dependency_type);

                auto *vector_FEEval_ptr =
                  vector_vars_map.at(dependency_index).at(dependency_type).get();
                vector_FEEval_ptr->reinit(cell);

                if (eval_flag_set.find(pair) != eval_flag_set.end())
                  {
                    const unsigned int &local_index = global_to_local_solution.at(pair);

                    Assert(src.size() > local_index,
                           dealii::ExcMessage(
                             "The provided src vector's size is below the given local "
                             "index = " +
                             std::to_string(local_index) +
                             " for global index = " + std::to_string(dependency_index) +
                             "  and type = " + to_string(dependency_type)));

                    vector_FEEval_ptr->read_dof_values_plain(*(src.at(local_index)));
                    vector_FEEval_ptr->evaluate(eval_flag_set.at(pair));
                  }
              }
          }
      }
  };

  if (solve_type == solveType::EXPLICIT_RHS || solve_type == solveType::POSTPROCESS)
    {
      reinit_and_eval_map(subset_attributes.begin()->second.eval_flag_set_RHS,
                          subset_attributes.begin()->second.dependency_set_RHS);
      return;
    }
  if (src.empty())
    {
      return;
    }

  Assert(subset_attributes.size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      reinit_and_eval_map(subset_attributes.begin()->second.eval_flag_set_LHS,
                          subset_attributes.begin()->second.dependency_set_LHS);
    }
  else
    {
      reinit_and_eval_map(subset_attributes.begin()->second.eval_flag_set_RHS,
                          subset_attributes.begin()->second.dependency_set_RHS);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::reinit_and_eval(const VectorType &src,
                                                        unsigned int      cell)
{
  // Reinit and eval values for the given dependency set. Note the dependency set includes
  // the variable we're evaluating, which may or may not be an actually dependency. For
  // this reason, I selectively read dofs and evaluate the flags.
  auto reinit_and_eval_map =
    [&](const std::unordered_map<std::pair<unsigned int, dependencyType>,
                                 dealii::EvaluationFlags::EvaluationFlags,
                                 pairHash>                                &eval_flag_set,
        const std::map<unsigned int, std::map<dependencyType, fieldType>> &dependency_set)
  {
    for (const auto &[dependency_index, map] : dependency_set)
      {
        for (const auto &[dependency_type, field_type] : map)
          {
            const auto &pair = std::make_pair(dependency_index, dependency_type);

            Assert(global_to_local_solution.find(pair) != global_to_local_solution.end(),
                   dealii::ExcMessage(
                     "The global to local mapping does not exists for global index = " +
                     std::to_string(dependency_index) +
                     "  and type = " + to_string(dependency_type)));

            if (field_type == fieldType::SCALAR)
              {
                scalar_FEEval_exists(dependency_index, dependency_type);

                auto *scalar_FEEval_ptr =
                  scalar_vars_map.at(dependency_index).at(dependency_type).get();
                scalar_FEEval_ptr->reinit(cell);
                scalar_FEEval_ptr->read_dof_values_plain(src);
                scalar_FEEval_ptr->evaluate(eval_flag_set.at(pair));
              }
            else
              {
                vector_FEEval_exists(dependency_index, dependency_type);

                auto *vector_FEEval_ptr =
                  vector_vars_map.at(dependency_index).at(dependency_type).get();
                vector_FEEval_ptr->reinit(cell);
                vector_FEEval_ptr->read_dof_values_plain(src);
                vector_FEEval_ptr->evaluate(eval_flag_set.at(pair));
              }
          }
      }
  };

  Assert(subset_attributes.size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      reinit_and_eval_map(subset_attributes.begin()->second.eval_flag_set_LHS,
                          subset_attributes.begin()->second.dependency_set_LHS);
    }
  else
    {
      Assert(false,
             dealii::ExcMessage(
               "reinit_and_eval(src) should only be called for LHS evaluations"));
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::reinit(unsigned int        cell,
                                               const unsigned int &global_variable_index)
{
  auto reinit_map = [&](const std::unordered_map<std::pair<unsigned int, dependencyType>,
                                                 dealii::EvaluationFlags::EvaluationFlags,
                                                 pairHash> &eval_flag_set)
  {
    for (const auto &[pair, flags] : eval_flag_set)
      {
        const unsigned int   &dependency_index = pair.first;
        const dependencyType &dependency_type  = pair.second;

        Assert(subset_attributes.find(dependency_index) != subset_attributes.end(),
               dealii::ExcMessage(
                 "The subset attribute entry does not exists for global index = " +
                 std::to_string(dependency_index)));

        if (subset_attributes.at(dependency_index).field_type == fieldType::SCALAR)
          {
            scalar_FEEval_exists(dependency_index, dependency_type);

            auto *scalar_FEEval_ptr =
              scalar_vars_map.at(dependency_index).at(dependency_type).get();
            scalar_FEEval_ptr->reinit(cell);
          }
        else
          {
            vector_FEEval_exists(dependency_index, dependency_type);

            auto *vector_FEEval_ptr =
              vector_vars_map.at(dependency_index).at(dependency_type).get();
            vector_FEEval_ptr->reinit(cell);
          }
      }
  };

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      Assert(subset_attributes.find(global_variable_index) != subset_attributes.end(),
             dealii::ExcMessage(
               "The subset attribute entry does not exists for global index = " +
               std::to_string(global_variable_index)));

      reinit_map(subset_attributes.at(global_variable_index).eval_flag_set_LHS);
    }
  else
    {
      Assert(subset_attributes.find(global_variable_index) != subset_attributes.end(),
             dealii::ExcMessage(
               "The subset attribute entry does not exists for global index = " +
               std::to_string(global_variable_index)));

      reinit_map(subset_attributes.at(global_variable_index).eval_flag_set_RHS);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval(const unsigned int &global_variable_index)
{
  auto eval_map = [&](const std::unordered_map<std::pair<unsigned int, dependencyType>,
                                               dealii::EvaluationFlags::EvaluationFlags,
                                               pairHash> &eval_flag_set)
  {
    for (const auto &[pair, flags] : eval_flag_set)
      {
        const unsigned int   &dependency_index = pair.first;
        const dependencyType &dependency_type  = pair.second;

        Assert(subset_attributes.find(dependency_index) != subset_attributes.end(),
               dealii::ExcMessage(
                 "The subset attribute entry does not exists for global index = " +
                 std::to_string(dependency_index)));

        if (subset_attributes.at(dependency_index).field_type == fieldType::SCALAR)
          {
            scalar_FEEval_exists(dependency_index, dependency_type);

            auto *scalar_FEEval_ptr =
              scalar_vars_map.at(dependency_index).at(dependency_type).get();
            scalar_FEEval_ptr->evaluate(flags);
          }
        else
          {
            vector_FEEval_exists(dependency_index, dependency_type);

            auto *vector_FEEval_ptr =
              vector_vars_map.at(dependency_index).at(dependency_type).get();
            vector_FEEval_ptr->evaluate(flags);
          }
      }
  };

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      Assert(subset_attributes.find(global_variable_index) != subset_attributes.end(),
             dealii::ExcMessage(
               "The subset attribute entry does not exists for global index = " +
               std::to_string(global_variable_index)));

      eval_map(subset_attributes.at(global_variable_index).eval_flag_set_LHS);
    }
  else
    {
      Assert(subset_attributes.find(global_variable_index) != subset_attributes.end(),
             dealii::ExcMessage(
               "The subset attribute entry does not exists for global index = " +
               std::to_string(global_variable_index)));

      eval_map(subset_attributes.at(global_variable_index).eval_flag_set_RHS);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::integrate(
  const unsigned int &global_variable_index)
{
  Assert(subset_attributes.find(global_variable_index) != subset_attributes.end(),
         dealii::ExcMessage(
           "The subset attribute entry does not exists for global index = " +
           std::to_string(global_variable_index)));

  const auto &variable = subset_attributes.at(global_variable_index);

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      if (variable.field_type == fieldType::SCALAR)
        {
          scalar_FEEval_exists(global_variable_index, dependencyType::CHANGE);

          auto *scalar_FEEval_ptr =
            scalar_vars_map.at(global_variable_index).at(dependencyType::CHANGE).get();
          scalar_FEEval_ptr->integrate(variable.eval_flags_residual_LHS);
        }
      else
        {
          vector_FEEval_exists(global_variable_index, dependencyType::CHANGE);

          auto *vector_FEEval_ptr =
            vector_vars_map.at(global_variable_index).at(dependencyType::CHANGE).get();
          vector_FEEval_ptr->integrate(variable.eval_flags_residual_LHS);
        }
    }
  else
    {
      AssertThrow(false,
                  dealii::ExcMessage(
                    "Integrate called for a solve type that is not NONEXPLICIT_LHS."));
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::integrate_and_distribute(
  std::vector<VectorType *> &dst)
{
  auto integrate_and_distribute_map =
    [&](const dealii::EvaluationFlags::EvaluationFlags &residual_flag_set,
        const dependencyType                           &dependency_type,
        const unsigned int                             &residual_index)
  {
    Assert(
      global_to_local_solution.find(std::make_pair(residual_index, dependency_type)) !=
        global_to_local_solution.end(),
      dealii::ExcMessage(
        "The global to local mapping does not exists for global index = " +
        std::to_string(residual_index) + "  and type = " + to_string(dependency_type)));

    const unsigned int &local_index =
      global_to_local_solution.at(std::make_pair(residual_index, dependency_type));

    Assert(dst.size() > local_index,
           dealii::ExcMessage(
             "The provided dst vector's size is below the given local index = " +
             std::to_string(local_index) +
             " for global index = " + std::to_string(residual_index) +
             "  and type = " + to_string(dependency_type)));
    Assert(subset_attributes.find(residual_index) != subset_attributes.end(),
           dealii::ExcMessage(
             "The subset attribute entry does not exists for global index = " +
             std::to_string(residual_index)));

    if (subset_attributes.at(residual_index).field_type == fieldType::SCALAR)
      {
        scalar_FEEval_exists(residual_index, dependency_type);

        auto *scalar_FEEval_ptr =
          scalar_vars_map.at(residual_index).at(dependency_type).get();
        scalar_FEEval_ptr->integrate_scatter(residual_flag_set, *(dst.at(local_index)));
      }
    else
      {
        vector_FEEval_exists(residual_index, dependency_type);

        auto *vector_FEEval_ptr =
          vector_vars_map.at(residual_index).at(dependency_type).get();
        vector_FEEval_ptr->integrate_scatter(residual_flag_set, *(dst.at(local_index)));
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
void
variableContainer<dim, degree, number>::integrate_and_distribute(VectorType &dst)
{
  auto integrate_and_distribute_map =
    [&](const dealii::EvaluationFlags::EvaluationFlags &residual_flag_set,
        const dependencyType                           &dependency_type,
        const unsigned int                             &residual_index)
  {
    Assert(subset_attributes.find(residual_index) != subset_attributes.end(),
           dealii::ExcMessage(
             "The subset attribute entry does not exists for global index = " +
             std::to_string(residual_index)));

    if (subset_attributes.at(residual_index).field_type == fieldType::SCALAR)
      {
        scalar_FEEval_exists(residual_index, dependency_type);

        auto *scalar_FEEval_ptr =
          scalar_vars_map.at(residual_index).at(dependency_type).get();
        scalar_FEEval_ptr->integrate_scatter(residual_flag_set, dst);
      }
    else
      {
        vector_FEEval_exists(residual_index, dependency_type);

        auto *vector_FEEval_ptr =
          vector_vars_map.at(residual_index).at(dependency_type).get();
        vector_FEEval_ptr->integrate_scatter(residual_flag_set, dst);
      }
  };

  Assert(subset_attributes.size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      integrate_and_distribute_map(
        subset_attributes.begin()->second.eval_flags_residual_LHS,
        dependencyType::CHANGE,
        subset_attributes.begin()->first);
    }
  else
    {
      integrate_and_distribute_map(
        subset_attributes.begin()->second.eval_flags_residual_RHS,
        dependencyType::NORMAL,
        subset_attributes.begin()->first);
    }
}

template <int dim, int degree, typename number>
typename variableContainer<dim, degree, number>::size_type
variableContainer<dim, degree, number>::get_scalar_value(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::values);
  scalar_FEEval_exists(global_variable_index, dependency_type);
#endif

  return scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_value(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_scalar_gradient(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  scalar_FEEval_exists(global_variable_index, dependency_type);
#endif

  return scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_gradient(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_scalar_hessian(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  scalar_FEEval_exists(global_variable_index, dependency_type);
#endif

  return scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_hessian(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_scalar_hessian_diagonal(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  scalar_FEEval_exists(global_variable_index, dependency_type);
#endif

  return scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_hessian_diagonal(q_point);
}

template <int dim, int degree, typename number>
typename variableContainer<dim, degree, number>::size_type
variableContainer<dim, degree, number>::get_scalar_laplacian(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  scalar_FEEval_exists(global_variable_index, dependency_type);
#endif

  return scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_laplacian(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_value(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::values);
  vector_FEEval_exists(global_variable_index, dependency_type);
#endif

  const auto &value =
    vector_vars_map.at(global_variable_index).at(dependency_type)->get_value(q_point);

  if constexpr (dim == 1)
    {
      // Wrap the value for consistency
      dealii::Tensor<1, dim, size_type> wrapper;
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
dealii::Tensor<2, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_gradient(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  vector_FEEval_exists(global_variable_index, dependency_type);
#endif

  const auto &grad =
    vector_vars_map.at(global_variable_index).at(dependency_type)->get_gradient(q_point);

  if constexpr (dim == 1)
    {
      // Wrap the value for consistency
      dealii::Tensor<2, dim, size_type> wrapper;
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
dealii::Tensor<3, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_hessian(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  vector_FEEval_exists(global_variable_index, dependency_type);
#endif

  const auto &hess =
    vector_vars_map.at(global_variable_index).at(dependency_type)->get_hessian(q_point);

  if constexpr (dim == 1)
    {
      // Wrap the value for consistency
      dealii::Tensor<3, dim, size_type> wrapper;
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
dealii::Tensor<2, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_hessian_diagonal(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  vector_FEEval_exists(global_variable_index, dependency_type);
#endif

  const auto &hess_diag = vector_vars_map.at(global_variable_index)
                            .at(dependency_type)
                            ->get_hessian_diagonal(q_point);

  if constexpr (dim == 1)
    {
      // Wrap the value for consistency
      dealii::Tensor<2, dim, size_type> wrapper;
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
dealii::Tensor<1, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_laplacian(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  vector_FEEval_exists(global_variable_index, dependency_type);
#endif

  const auto &lap =
    vector_vars_map.at(global_variable_index).at(dependency_type)->get_laplacian(q_point);

  if constexpr (dim == 1)
    {
      // Wrap the value for consistency
      dealii::Tensor<1, dim, size_type> wrapper;
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
typename variableContainer<dim, degree, number>::size_type
variableContainer<dim, degree, number>::get_vector_divergence(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  vector_FEEval_exists(global_variable_index, dependency_type);
#endif

  return vector_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_divergence(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_symmetric_gradient(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  vector_FEEval_exists(global_variable_index, dependency_type);
#endif

  return vector_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->get_symmetric_gradient(q_point);
}

template <int dim, int degree, typename number>
dealii::Tensor<1,
               (dim == 2 ? 1 : dim),
               typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_curl(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  vector_FEEval_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      Assert(false, dealii::ExcMessage("PRISMS-PF Error: Curl is nonsensical for 1D."));
      return dealii::Tensor<1, (dim == 2 ? 1 : dim), size_type> {};
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
  const unsigned int   &global_variable_index,
  const size_type      &val,
  const dependencyType &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  scalar_FEEval_exists(global_variable_index, dependency_type);
#endif

  scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->submit_value(val, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_scalar_gradient_term(
  const unsigned int                      &global_variable_index,
  const dealii::Tensor<1, dim, size_type> &grad,
  const dependencyType                    &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  scalar_FEEval_exists(global_variable_index, dependency_type);
#endif

  scalar_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->submit_gradient(grad, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_vector_value_term(
  const unsigned int                      &global_variable_index,
  const dealii::Tensor<1, dim, size_type> &val,
  const dependencyType                    &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  vector_FEEval_exists(global_variable_index, dependency_type);

#endif

  vector_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->submit_value(val, q_point);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_vector_gradient_term(
  const unsigned int                      &global_variable_index,
  const dealii::Tensor<2, dim, size_type> &grad,
  const dependencyType                    &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  vector_FEEval_exists(global_variable_index, dependency_type);
#endif

  vector_vars_map.at(global_variable_index)
    .at(dependency_type)
    ->submit_gradient(grad, q_point);
}

INSTANTIATE_TRI_TEMPLATE(variableContainer)

PRISMS_PF_END_NAMESPACE