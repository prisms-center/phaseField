// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <map>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class handles the generation and representation of each field
 * with an associated index.
 *
 * As the complexity of the types of fields we can have grows, it becomes
 * burdensome to keep track of nest maps and/or maps with pairs and tuples. This
 * class aims to unify the mapping of unique fields to some index so that we can
 * use vector to manipulate data. Additionally, it centralizes the index
 * mapping, as opposed to having local mappings spread throughout the code.
 */

class indexMap
{
public:
  using map_index = unsigned int;

  /**
   *  \brief Blank constructor.
   */
  indexMap() = default;

  /**
   * \brief Default Constructor.
   */
  explicit indexMap(
    const std::map<unsigned int, variableAttributes> &_variable_attributes);

  /**
   * \brief Initialize the mapping.
   */
  void
  init();

private:
  /**
   * \brief Variable attributes.
   */
  const std::map<unsigned int, variableAttributes> *variable_attributes = nullptr;

  /**
   * \brief Mapping between field index, dependency type, and the numbering in
   * indexMap.
   */
  std::vector<std::pair<unsigned int, dependencyType>> mapping;
};

indexMap::indexMap(const std::map<unsigned int, variableAttributes> &_variable_attributes)
  : variable_attributes(&_variable_attributes)
{
  init();
}

inline void
indexMap::init()
{
  Assert(variable_attributes != nullptr, dealii::ExcNotInitialized());

  // Loop through the variable_attributes and fill in the mappings. This is done
  // in such a way that the vector is ordered in terms of index and then
  // dependencyType. For example, 1 NORMAL, 1 OLD_1, 2 NORMAL, 3 NORMAL, 3
  // OLD_1, 3 OLD_2...
}

PRISMS_PF_END_NAMESPACE
