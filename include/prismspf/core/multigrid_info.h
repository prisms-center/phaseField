// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/linear_solve_parameters.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <map>
#include <set>
#include <tuple>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class centralized the collection of multigrid data, which is of particular
 * important when dealing with multiple fields that have multigrid enabled.
 */
template <unsigned int dim>
class MGInfo
{
public:
  using Min = unsigned int;
  using Max = unsigned int;

  /**
   * \brief Constructor.
   */
  explicit MGInfo(const UserInputParameters<dim> &_user_inputs);

  /**
   * \brief If multigrid is enabled.
   */
  [[nodiscard]] bool
  has_multigrid() const
  {
    return multigrid_on;
  }

  /**
   * \brief Get the minimum multigrid level.
   */
  [[nodiscard]] Min
  get_mg_min_level() const
  {
    Assert(multigrid_on, dealii::ExcNotInitialized());
    return global_mg_level.first;
  }

  /**
   * \brief Get the maximum multigrid level.
   */
  [[nodiscard]] Max
  get_mg_max_level() const
  {
    Assert(multigrid_on, dealii::ExcNotInitialized());
    return global_mg_level.second;
  }

  /**
   * \brief Get the collection of minimum multigrid levels for the LHS fields.
   */
  [[nodiscard]] std::set<std::tuple<Types::Index, DependencyType, Min>>
  get_lhs_fields() const
  {
    Assert(multigrid_on, dealii::ExcNotInitialized());
    return lhs_fields;
  }

  /**
   * \brief Get the depth of the multigrid hierarchy
   */
  [[nodiscard]] unsigned int
  get_mg_depth() const
  {
    Assert(multigrid_on, dealii::ExcNotInitialized());
    return mg_levels.size();
  }

  /**
   * \brief Get the depth of the multigrid hierarchy
   */
  [[nodiscard]] unsigned int
  get_mg_breadth(unsigned int level) const
  {
    Assert(multigrid_on, dealii::ExcNotInitialized());
    Assert(level < mg_levels.size(),
           dealii::ExcMessage(
             "The multigrid hierarchy does not contain the relative level = " +
             std::to_string(level)));
    return mg_levels[level].size();
  }

  /**
   * \brief Get the global index for a given local index and relative level.
   */
  [[nodiscard]] unsigned int
  get_global_index(unsigned int index, unsigned int level) const
  {
    Assert(multigrid_on, dealii::ExcNotInitialized());
    Assert(level < mg_levels.size(),
           dealii::ExcMessage(
             "The multigrid hierarchy does not contain the relative level = " +
             std::to_string(level)));
    Assert(index < mg_levels[level].size(),
           dealii::ExcMessage(
             "The multigrid hierarchy does not contain the local index = " +
             std::to_string(index)));
    return mg_levels[level][index];
  }

  /**
   * \brief Get the local index for a given global index and level.
   */
  [[nodiscard]] unsigned int
  get_local_index(unsigned int global_index, unsigned int level) const
  {
    Assert(multigrid_on, dealii::ExcNotInitialized());
    const unsigned int relative_level = level - global_mg_level.first;
    Assert(relative_level < mg_levels.size(),
           dealii::ExcMessage("The multigrid hierarchy does not contain level = " +
                              std::to_string(level)));
    auto iterator = std::find(mg_levels[relative_level].begin(),
                              mg_levels[relative_level].end(),
                              global_index);
    Assert(iterator != mg_levels[relative_level].end(),
           dealii::ExcMessage("The multigrid level does not contain index = " +
                              std::to_string(global_index)));
    return std::distance(mg_levels[relative_level].begin(), iterator);
  }

  /**
   * \brief Print information about the instances of this class.
   */
  void
  print();

private:
  /**
   * \brief User inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * \brief Whether multigrid is enabled.
   */
  bool multigrid_on = false;

  /**
   * \brief Global Min and Max multigrid levels.
   */
  std::pair<Min, Max> global_mg_level;

  /**
   * \brief The collection of LHS fields and their minimum multigrid level that need to be
   * initialized for the LHS of gmg fields.
   */
  std::set<std::tuple<Types::Index, DependencyType, Min>> lhs_fields;

  /**
   * \brief Vector hierarchy of required multigrid levels.
   */
  std::vector<std::vector<unsigned int>> mg_levels;
};

template <unsigned int dim>
inline MGInfo<dim>::MGInfo(const UserInputParameters<dim> &_user_inputs)
  : user_inputs(&_user_inputs)
{
  const Max max_mg_level =
    user_inputs->get_spatial_discretization().get_has_adaptivity()
      ? user_inputs->get_spatial_discretization().get_max_refinement()
      : user_inputs->get_spatial_discretization().get_global_refinement();
  Min min_mg_level = UINT_MAX;

  std::set<Types::Index>      fields_with_multigrid;
  std::map<Types::Index, Min> min_levels;
  for (const auto &[index, linear_solver] :
       user_inputs->get_linear_solve_parameters().get_linear_solve_parameters())
    {
      if (linear_solver.preconditioner == PreconditionerType::GMG)
        {
          const Min min_level = linear_solver.min_mg_level;
          fields_with_multigrid.insert(index);
          min_levels.emplace(index, min_level);
          min_mg_level = std::min(min_mg_level, min_level);
        }
    }

  global_mg_level = std::make_pair(min_mg_level, max_mg_level);

  if (fields_with_multigrid.empty())
    {
      multigrid_on = false;
      return;
    }
  multigrid_on = true;

  // Create the set of tuples that let us know which fields have what multigrid levels.
  // Note that we don't have to do this if the DependencyType is Change, although we must
  // always include DependencyType::Normal.
  std::set<std::tuple<Types::Index, DependencyType, Min>> all_lhs_fields;
  for (const auto &index : fields_with_multigrid)
    {
      const auto &variable = user_inputs->get_variable_attributes().at(index);

      // First add the DependencyType::Change
      all_lhs_fields.insert(
        std::make_tuple(index, DependencyType::Change, min_levels.at(index)));

      // Then add the LHS dependencies. For now I just add all of them and go back to trim
      // the duplicates that have the same variable but different Min multigrid levels.
      for (const auto &[pair, eval_flag] : variable.get_eval_flag_set_lhs())
        {
          // Skip if the eval flags is not set (e.i., nothing)
          if (eval_flag == 0U)
            {
              continue;
            }

          AssertThrow(pair.second == DependencyType::Normal ||
                        pair.second == DependencyType::Change,
                      dealii::ExcNotImplemented());

          all_lhs_fields.insert(
            std::make_tuple(pair.first, pair.second, min_levels.at(index)));
        }
    }
  // Trim fields that have the same first and second entry of the tuple.
  std::map<std::pair<Types::Index, DependencyType>, Min> lowest_mg_level;
  for (const auto &tuple : all_lhs_fields)
    {
      auto      key       = std::make_pair(std::get<0>(tuple), std::get<1>(tuple));
      const Min min_level = std::get<2>(tuple);

      auto iterator = lowest_mg_level.find(key);
      if (iterator == lowest_mg_level.end() || min_level < iterator->second)
        {
          lowest_mg_level[key] = min_level;
        }
    }
  for (const auto &[key, min_level] : lowest_mg_level)
    {
      lhs_fields.insert(std::make_tuple(key.first, key.second, min_level));
    }

  mg_levels.resize(global_mg_level.second - global_mg_level.first + 1);
  for (const auto &[index, dependency, min_level] : lhs_fields)
    {
      const unsigned int relative_level = min_level - global_mg_level.first;
      for (unsigned int level = relative_level; level < mg_levels.size(); level++)
        {
          mg_levels[level].push_back(index);
        }
    }
}

template <unsigned int dim>
void
MGInfo<dim>::print()
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
    << "  MGInfo\n"
    << "================================================\n";
  if (!multigrid_on)
    {
      ConditionalOStreams::pout_summary()
        << "  There are no fields with multigrid enabled\n\n"
        << std::flush;
      return;
    }
  ConditionalOStreams::pout_summary() << "  Global min = " << global_mg_level.first
                                      << " and max = " << global_mg_level.second << "\n"
                                      << "  LHS dependency fields:\n";
  for (const auto &[field, dependency, minimum] : lhs_fields)
    {
      ConditionalOStreams::pout_summary()
        << "    Index " << field << " Type " << to_string(dependency) << " Min "
        << minimum << "\n";
    }
  unsigned int outer_index = 0;
  ConditionalOStreams::pout_summary() << "  LHS index hierarchy:\n";
  for (const auto &vector : mg_levels)
    {
      ConditionalOStreams::pout_summary() << "    Outer index " << outer_index << "\n";
      unsigned int inner_index = 0;
      for (const auto &index : vector)
        {
          ConditionalOStreams::pout_summary()
            << "      Inner index " << inner_index << " Field index: " << index << "\n";
          inner_index++;
        }
      outer_index++;
    }
  ConditionalOStreams::pout_summary() << "\n" << std::flush;
}

PRISMS_PF_END_NAMESPACE
