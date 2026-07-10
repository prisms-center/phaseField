#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>

#include <prismspf/core/exceptions.h>

#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/config.h>

#include <cfloat>

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

  std::vector<std::string_view> subsection_stack;
  std::vector<Parameter>        parameters;
};

PRISMS_PF_END_NAMESPACE
