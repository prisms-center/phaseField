#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

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
      parameter_handler.declare_alias(existing_entry_name, alias);
      declare_aliases(parameter_handler, existing_entry_name, generate_aliases(alias));
    }
}

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
