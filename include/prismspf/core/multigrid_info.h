// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

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
template <int dim>
class MGInfo
{
public:
  using min = unsigned int;
  using max = unsigned int;

  /**
   * \brief Constructor.
   */
  explicit MGInfo(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Print information about the instances of this class.
   */
  void
  print();

private:
  /**
   * \brief User inputs.
   */
  const userInputParameters<dim> *user_inputs;

  /**
   * \brief Whether multigrid is enabled.
   */
  bool has_multigrid = false;

  /**
   * \brief Global min and max multigrid levels.
   */
  std::pair<min, max> global_mg_level;

  /**
   * \brief The collection of LHS fields and their minimum multigrid level that need to be
   * initialized for the LHS of gmg fields.
   */
  std::set<std::tuple<types::index, dependencyType, min>> lhs_fields;
};

template <int dim>
inline MGInfo<dim>::MGInfo(const userInputParameters<dim> &_user_inputs)
  : user_inputs(&_user_inputs)
{
  const max max_mg_level = user_inputs->spatial_discretization.has_adaptivity
                             ? user_inputs->spatial_discretization.max_refinement
                             : user_inputs->spatial_discretization.global_refinement;
  min       min_mg_level = UINT_MAX;

  std::set<types::index>      fields_with_multigrid;
  std::map<types::index, min> min_levels;
  for (const auto &[index, linear_solver] :
       user_inputs->linear_solve_parameters.linear_solve)
    {
      if (linear_solver.preconditioner == preconditionerType::GMG)
        {
          const min min_level = linear_solver.min_mg_level;
          fields_with_multigrid.insert(index);
          min_levels.emplace(index, min_level);
          min_mg_level = std::min(min_mg_level, min_level);
        }
    }

  global_mg_level = std::make_pair(min_mg_level, max_mg_level);

  if (fields_with_multigrid.empty())
    {
      has_multigrid = false;
      return;
    }
  has_multigrid = true;

  // Create the set of tuples that let us know which fields have what multigrid levels.
  // Note that we don't have to do this if the dependencyType is CHANGE, although we must
  // always include dependencyType::NORMAL.
  std::set<std::tuple<types::index, dependencyType, min>> all_lhs_fields;
  for (const auto &index : fields_with_multigrid)
    {
      const auto &variable = user_inputs->var_attributes->at(index);

      // First add the dependencyType::CHANGE
      all_lhs_fields.insert(
        std::make_tuple(index, dependencyType::CHANGE, min_levels.at(index)));

      // Then add the LHS dependencies. For now I just add all of them and go back to trim
      // the duplicates that have the same variable but different min multigrid levels.
      for (const auto &[pair, eval_flag] : variable.eval_flag_set_LHS)
        {
          // Skip if the eval flags is not set (e.i., nothing)
          if (eval_flag == 0U)
            {
              continue;
            }

          all_lhs_fields.insert(
            std::make_tuple(pair.first, pair.second, min_levels.at(index)));
        }
    }
  // Trim fields that have the same first and second entry of the tuple.
  std::map<std::pair<types::index, dependencyType>, min> lowest_mg_level;
  for (const auto &tuple : all_lhs_fields)
    {
      auto key       = std::make_pair(std::get<0>(tuple), std::get<1>(tuple));
      min  min_level = std::get<2>(tuple);

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
}

template <int dim>
void
MGInfo<dim>::print()
{
  conditionalOStreams::pout_base()
    << "================================================\n"
    << "  MGInfo\n"
    << "================================================\n";
  if (!has_multigrid)
    {
      conditionalOStreams::pout_base()
        << "  There are no fields with multigrid enabled\n\n"
        << std::flush;
      return;
    }
  conditionalOStreams::pout_base() << "  Global min = " << global_mg_level.first
                                   << " and max = " << global_mg_level.second << "\n"
                                   << "  LHS dependency fields:\n";
  for (const auto &[field, dependency, minimum] : lhs_fields)
    {
      conditionalOStreams::pout_base()
        << "    Index " << field << " Type " << to_string(dependency) << " Min "
        << minimum << "\n";
    }
  conditionalOStreams::pout_base() << std::flush;
}

PRISMS_PF_END_NAMESPACE
