#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/utilities.h>

#include <prismspf/core/exceptions.h>

#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/config.h>

#include <cfloat>
#include <concepts>
#include <ranges>

PRISMS_PF_BEGIN_NAMESPACE

enum class Pattern : std::uint8_t
{
  Anything,
  Selection,
  Bool,
  Integer,
  PositiveInteger,
  NegativeInteger,
  Double,
  PositiveDouble,
  NegativeDouble,
  StringList,
  PositiveIntegerList
};

std::string_view
to_string(Pattern pattern)
{
  switch (pattern)
    {
      case Pattern::Anything:
        return "Anything";
      case Pattern::Selection:
        return "Selection";
      case Pattern::Bool:
        return "Bool";
      case Pattern::Integer:
        return "Integer";
      case Pattern::PositiveInteger:
        return "PositiveInteger";
      case Pattern::NegativeInteger:
        return "NegativeInteger";
      case Pattern::Double:
        return "Double";
      case Pattern::PositiveDouble:
        return "PositiveDouble";
      case Pattern::NegativeDouble:
        return "NegativeDouble";
      case Pattern::StringList:
        return "StringList";
      case Pattern::PositiveIntegerList:
        return "PositiveIntegerList";
    }

  return "Unknown";
}

struct ParameterData
{
  using String = std::string_view;

  String            entry;
  String            default_value;
  Pattern           pattern;
  String            pattern_options;
  String            documentation;
  String            warnings;
  String            category;
  String            example;
  std::list<String> aliases;
  bool              advanced   = false;
  bool              deprecated = false;
};

std::string
to_html_row(const ParameterData &p)
{
  std::ostringstream html;

  html << "<tr>";
  html << "<td>" << p.entry << "</td>";
  html << "<td>" << p.default_value << "</td>";
  html << "<td>" << to_string(p.pattern) << "</td>";
  html << "<td>" << p.pattern_options << "</td>";
  html << "<td>" << p.documentation << "</td>";
  html << "<td>" << p.warnings << "</td>";
  html << "<td>" << p.category << "</td>";
  html << "<td>" << p.example << "</td>";

  html << "<td>";
  for (const auto &alias : p.aliases)
    {
      html << alias << "<br>";
    }
  html << "</td>";

  html << "<td>" << (p.advanced ? "true" : "false") << "</td>";
  html << "<td>" << (p.deprecated ? "true" : "false") << "</td>";

  html << "</tr>";

  return html.str();
}

struct Parameter
{
  std::vector<std::string_view> subsection_path;
  ParameterData                 data;
};

namespace Patterns
{
  // NOLINTBEGIN(readability-identifier-naming)

  inline dealii::Patterns::Anything
  Anything()
  {
    return {};
  }

  inline dealii::Patterns::Selection
  Selection(std::string_view options)
  {
    return {std::string(options)};
  }

  inline dealii::Patterns::Bool
  Bool()
  {
    return {};
  }

  inline dealii::Patterns::Integer
  Integer()
  {
    return {-INT_MAX, INT_MAX};
  }

  inline dealii::Patterns::Integer
  PositiveInteger()
  {
    return {0, INT_MAX};
  }

  inline dealii::Patterns::Integer
  NegativeInteger()
  {
    return {-INT_MAX, 0};
  }

  inline dealii::Patterns::Double
  Double()
  {
    return {-DBL_MAX, DBL_MAX};
  }

  inline dealii::Patterns::Double
  PositiveDouble()
  {
    return {0.0, DBL_MAX};
  }

  inline dealii::Patterns::Double
  NegativeDouble()
  {
    return {-DBL_MAX, 0.0};
  }

  inline dealii::Patterns::List
  StringList()
  {
    return {Anything(), 0, INT_MAX, ","};
  }

  inline dealii::Patterns::List
  PositiveIntegerList()
  {
    return {PositiveInteger(), 0, INT_MAX, ","};
  }

  // NOLINTEND(readability-identifier-naming)

} // namespace Patterns

/**
 * @brief Wrapper for deal.II's ParameterHandler class.
 *
 * This is a wrapper for deal.II's ParameterHandler class. We do this so we can provide a
 * little extra documentation and auto generate our parameter documentation from code.
 */
class ParameterHandler
{
  using String = std::string_view;

  void
  declare_entry(const ParameterData &data)
  {
    parameters.push_back({subsection_stack, data});

    auto declare = [&](const auto &pattern)
      {
        parameter_handler.declare_entry(std::string(data.entry),
                                        std::string(data.default_value),
                                        pattern);
      };

    switch (data.pattern)
      {
        case Pattern::Anything:
          declare(Patterns::Anything());
          break;
        case Pattern::Selection:
          declare(Patterns::Selection(data.pattern_options));
          break;
        case Pattern::Bool:
          declare(Patterns::Bool());
          break;
        case Pattern::Integer:
          declare(Patterns::Integer());
          break;
        case Pattern::PositiveInteger:
          declare(Patterns::PositiveInteger());
          break;
        case Pattern::NegativeInteger:
          declare(Patterns::NegativeInteger());
          break;
        case Pattern::Double:
          declare(Patterns::Double());
          break;
        case Pattern::PositiveDouble:
          declare(Patterns::PositiveDouble());
          break;
        case Pattern::NegativeDouble:
          declare(Patterns::NegativeDouble());
          break;
        case Pattern::StringList:
          declare(Patterns::StringList());
          break;
        case Pattern::PositiveIntegerList:
          declare(Patterns::PositiveIntegerList());
          break;
        default:
          AssertThrow(false, UnreachableCode());
      }

    declare_aliases_with_generated(parameter_handler,
                                   std::string(data.entry),
                                   data.aliases);
  }

  std::string
  get(const String &entry) const
  {
    return parameter_handler.get(std::string(entry));
  }

  int
  get_int(const String &entry) const
  {
    // TODO: Add check for safe casting
    return (int) parameter_handler.get_integer(std::string(entry));
  }

  unsigned int
  get_unsigned_int(const String &entry) const
  {
    // TODO: Add check for safe casting
    return (unsigned int) parameter_handler.get_integer(std::string(entry));
  }

  double
  get_double(const String &entry) const
  {
    return parameter_handler.get_double(std::string(entry));
  }

  bool
  get_bool(const String &entry) const
  {
    return parameter_handler.get_bool(std::string(entry));
  }

  template <typename Map>
  requires requires(Map &map, const typename Map::key_type &key) {
    map.find(key);
    map.end();
  }
  typename Map::mapped_type
  get_selection(const String &entry, Map &map) const
  {
    const auto value = get(entry);
    const auto iter  = map.find(value);
    AssertThrow(iter != map.end(), dealii::ExcMessage("Invalid selection"));
    return iter->second;
  }

  std::vector<std::string>
  get_string_list(const String &entry) const
  {
    return dealii::Utilities::split_string_list(
      parameter_handler.get(std::string(entry)));
  }

  std::vector<unsigned int>
  get_unsigned_int_list(const String &entry) const
  {
    // TODO: Add check for safe casting
    std::vector<unsigned int> list;

    for (const auto &list_entry : get_string_list(entry))
      {
        list.push_back((unsigned int) std::stoul(list_entry));
      }

    return list;
  }

  void
  enter_subsection(const String &subsection)
  {
    subsection_stack.emplace_back(subsection);
    parameter_handler.enter_subsection(std::string(subsection));
  }

  void
  leave_subsection()
  {
    subsection_stack.pop_back();
    parameter_handler.leave_subsection();
  }

  void
  clear_but_preserve_documentation()
  {
    parameter_handler.clear();
    subsection_stack.clear();
  }

  void
  clear()
  {
    parameter_handler.clear();
    subsection_stack.clear();
    parameters.clear();
  }

private:
  dealii::ParameterHandler parameter_handler;

  std::vector<String>    subsection_stack;
  std::vector<Parameter> parameters;
};

PRISMS_PF_END_NAMESPACE
