#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <cfloat>
#include <concepts>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

class SolveBlock;

/**
 * @brief Cartesian axis labels
 */
static constexpr std::array<std::string_view, 3> axis_labels {"x", "y", "z"};

/**
 * Declare multiple aliases for a parameter name.
 */
template <std::ranges::input_range Range>
requires std::convertible_to<std::ranges::range_reference_t<Range>, std::string_view>
static void
declare_aliases(dealii::ParameterHandler &parameter_handler,
                const std::string        &existing_entry_name,
                const Range              &aliases)
{
  for (const auto &alias : aliases)
    {
      parameter_handler.declare_alias(existing_entry_name, alias);
    }
}

/**
 * Generate common aliases for a parameter name.
 *
 * Given "user inputs", generates:
 *   - "user_inputs"
 *   - "User Inputs"
 *   - "User_Inputs"
 *   - "UserInputs"
 */
[[nodiscard]] static std::array<std::string, 4>
generate_aliases(std::string_view name)
{
  std::string snake;
  std::string title;
  std::string title_snake;
  std::string pascal;
  snake.reserve(name.size());
  title.reserve(name.size());
  title_snake.reserve(name.size());
  pascal.reserve(name.size());

  bool start_of_word = true;
  for (char character : name)
    {
      if (character == ' ')
        {
          snake += '_';
          title += ' ';
          title_snake += '_';
          start_of_word = true;
          continue;
        }

      const char lower = (char) std::tolower((unsigned char) character);
      const char upper = (char) std::toupper((unsigned char) character);

      snake += lower;
      title += start_of_word ? upper : lower;
      title_snake += start_of_word ? upper : lower;
      pascal += start_of_word ? upper : lower;

      start_of_word = false;
    }

  return {
    std::move(snake),
    std::move(title),
    std::move(title_snake),
    std::move(pascal),
  };
}

/**
 * Declare multiple aliases for a parameter name with generated aliases.
 */
template <std::ranges::input_range Range>
requires std::convertible_to<std::ranges::range_reference_t<Range>, std::string_view>
static void
declare_aliases_with_generated(dealii::ParameterHandler &parameter_handler,
                               const std::string        &existing_entry_name,
                               const Range              &aliases)
{
  declare_aliases(parameter_handler,
                  existing_entry_name,
                  generate_aliases(existing_entry_name));

  for (const auto &alias : aliases)
    {
      parameter_handler.declare_alias(existing_entry_name, std::string(alias));
      declare_aliases(parameter_handler, existing_entry_name, generate_aliases(alias));
    }
}

static void
declare_entry(dealii::ParameterHandler            &parameter_handler,
              const std::string_view              &entry,
              const std::string_view              &default_value,
              const dealii::Patterns::PatternBase &pattern,
              const std::string_view              &documentation,
              const std::list<std::string_view>   &aliases    = {},
              const std::string_view              &warnings   = {},
              const std::string_view              &example    = {},
              bool                                 deprecated = false)
{
  const auto _entry = std::string(entry);

  // Format the warning and example strings and append them to the documentation
  auto _doc = std::string(documentation);

  if (!warnings.empty())
    {
      _doc += "\n\n@warning " + std::string(warnings);
    }
  if (!example.empty())
    {
      _doc += "\n\n@example\n"
              "```text\n" +
              std::string(example) + "\n```\n";
    }

  parameter_handler.declare_entry(_entry,
                                  std::string(default_value),
                                  pattern,
                                  std::string(documentation));
  parameter_handler.mark_as_deprecated(_entry, deprecated);
  declare_aliases_with_generated(parameter_handler, _entry, aliases);
}

static std::string
get(const dealii::ParameterHandler &parameter_handler, const std::string_view &entry)
{
  return parameter_handler.get(std::string(entry));
}

static int
get_int(const dealii::ParameterHandler &parameter_handler, const std::string_view &entry)
{
  // TODO: Add check for safe casting
  return (int) parameter_handler.get_integer(std::string(entry));
}

static unsigned int
get_unsigned_int(const dealii::ParameterHandler &parameter_handler,
                 const std::string_view         &entry)
{
  // TODO: Add check for safe casting
  return (unsigned int) parameter_handler.get_integer(std::string(entry));
}

static double
get_double(const dealii::ParameterHandler &parameter_handler,
           const std::string_view         &entry)
{
  return parameter_handler.get_double(std::string(entry));
}

static bool
get_bool(const dealii::ParameterHandler &parameter_handler, const std::string_view &entry)
{
  return parameter_handler.get_bool(std::string(entry));
}

template <typename Map>
requires requires(Map &map, const typename Map::key_type &key) {
  map.find(key);
  map.end();
}
typename Map::mapped_type
get_selection(const dealii::ParameterHandler &parameter_handler,
              const std::string_view         &entry,
              Map                            &map)
{
  const auto value = get(parameter_handler, entry);
  const auto iter  = map.find(value);
  AssertThrow(iter != map.end(), dealii::ExcMessage("Invalid selection"));
  return iter->second;
}

static std::vector<std::string>
get_string_list(const dealii::ParameterHandler &parameter_handler,
                const std::string_view         &entry)
{
  return dealii::Utilities::split_string_list(parameter_handler.get(std::string(entry)));
}

static std::vector<unsigned int>
get_unsigned_int_list(const dealii::ParameterHandler &parameter_handler,
                      const std::string_view         &entry)
{
  // TODO: Add check for safe casting
  std::vector<unsigned int> list;

  for (const auto &list_entry : get_string_list(parameter_handler, entry))
    {
      list.push_back((unsigned int) std::stoul(list_entry));
    }

  return list;
}

namespace Patterns
{
  // NOLINTBEGIN(readability-identifier-naming)

  static inline dealii::Patterns::Anything
  Anything()
  {
    return {};
  }

  static inline dealii::Patterns::Selection
  Selection(std::string_view options)
  {
    return {std::string(options)};
  }

  static inline dealii::Patterns::Bool
  Bool()
  {
    return {};
  }

  static inline dealii::Patterns::Integer
  Integer()
  {
    return {-INT_MAX, INT_MAX};
  }

  static inline dealii::Patterns::Integer
  UnsignedInteger()
  {
    return {0, INT_MAX};
  }

  static inline dealii::Patterns::Integer
  PositiveInteger()
  {
    return {1, INT_MAX};
  }

  static inline dealii::Patterns::Integer
  NegativeInteger()
  {
    return {-INT_MAX, 0};
  }

  static inline dealii::Patterns::Double
  Double()
  {
    return {-DBL_MAX, DBL_MAX};
  }

  static inline dealii::Patterns::Double
  PositiveDouble()
  {
    return {0.0, DBL_MAX};
  }

  static inline dealii::Patterns::Double
  NegativeDouble()
  {
    return {-DBL_MAX, 0.0};
  }

  static inline dealii::Patterns::List
  StringList()
  {
    return {Anything(), 0, INT_MAX, ","};
  }

  static inline dealii::Patterns::List
  UnsignedIntegerList()
  {
    return {UnsignedInteger(), 0, INT_MAX, ","};
  }

  // NOLINTEND(readability-identifier-naming)

} // namespace Patterns

/**
 * @brief Virtual base class for parameter groups.
 *
 * This virtual base class ensures that we always have the following execution order.
 * 1. declare
 * 2. assign
 * 3. validate
 */
struct ParameterBase
{
  ParameterBase()          = default;
  virtual ~ParameterBase() = default;

  ParameterBase(const ParameterBase &) = default;
  ParameterBase(ParameterBase &&)      = default;

  ParameterBase &
  operator=(const ParameterBase &) = default;
  ParameterBase &
  operator=(ParameterBase &&) = default;

  /**
   * @brief Assign the parameters from file.
   */
  virtual void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              n_subsections = Numbers::default_subsections) = 0;

  /**
   * @brief Validate.
   */
  virtual void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const = 0;
};

PRISMS_PF_END_NAMESPACE
