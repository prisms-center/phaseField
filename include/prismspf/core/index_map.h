// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/config.h>

#include <map>
#include <type_traits>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <template <typename...> class T>
struct is_std_map : std::false_type
{};

template <>
struct is_std_map<std::map> : std::true_type
{};

template <typename T>
struct is_std_pair : std::false_type
{};

template <typename T, typename U>
struct is_std_pair<std::pair<T, U>> : std::true_type
{};

/**
 * \brief This class handlers the conversion of maps and nested maps into a linearly
 * indexed vector.
 *
 * For many places in the PRISMS-PF is makes sense to represented certain features with
 * maps. A great example is solution vectors. We have the field index and time the
 * solution field belongs too. This can either be repsented with a nested map or a map
 * between the vector and a pair. Of course, in places like variableContainer we end up
 * doing lots of lookups to the solution on every processor and over many points. This is
 * very computationally costly, albeit easy to read. For this reason, it is better to use
 * a vector.
 *
 * Passing some map to this class will construct some local indexing to use instead.
 */
template <template <typename...> class mapType,
          typename value,
          typename key1,
          typename key2 = key1>
class indexMap
{
public:
  using simpleMap = mapType<key1, value>;
  using nestedMap = mapType<key1, mapType<key2, value>>;

  /**
   * \brief Constructor.
   */
  explicit indexMap(const simpleMap &_map);

  /**
   * \brief Constructor with nested maps.
   */
  explicit indexMap(const nestedMap &_map);

private:
  /**
   * \brief Whether the mapType is a std::map.
   *
   * This is important because we want the same ordering regardless of whether a std::map
   * or std::unordered_map is passed. If a std::map is passed we don't have to spend time
   * creating a temp map to make sure the ordering is the same.
   */
  bool mapType_is_std_map = is_std_map<mapType>::value;

  /**
   * \brief Whether key1 is a pair.
   */
  bool key1_is_pair = is_std_pair<key1>::value;

  /**
   * \brief Whether key2 is a pair. This is unsupported and will immediately throw an
   * error.
   */
  bool key2_is_pair = is_std_pair<key2>::value;

  /**
   * \brief Vector of the keys.
   */
  std::vector<key1> vector_keys;

  /**
   * \brief Vector of the nested keys. This is only used for nested maps.
   */
  std::vector<std::pair<key1, key2>> vector_nested_keys;

  /**
   * \brief Vector of the values.
   */
  std::vector<const value *> vector_values;
};

PRISMS_PF_END_NAMESPACE