#pragma once

#include <prismspf/config.h>

#include <concepts>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <random>
#include <span>
#include <string>
#include <string_view>

PRISMS_PF_BEGIN_NAMESPACE

template <typename T>
concept StringLike = std::convertible_to<T, std::string_view>;

class TerminalColor
{
public:
  TerminalColor() = delete;

  // NOLINTBEGIN(readability-identifier-naming)

  static constexpr std::string_view RESET = "\033[0m";

  static constexpr std::string_view ERASE_SCREEN = "\033[2J";

  static constexpr std::string_view BOLD      = "\033[1m";
  static constexpr std::string_view DIM       = "\033[2m";
  static constexpr std::string_view ITALIC    = "\033[3m";
  static constexpr std::string_view UNDERLINE = "\033[4m";
  static constexpr std::string_view BLINK     = "\033[5m";
  static constexpr std::string_view STRIKE    = "\033[9m";

  static constexpr std::string_view BLACK   = "\033[30m";
  static constexpr std::string_view RED     = "\033[31m";
  static constexpr std::string_view GREEN   = "\033[32m";
  static constexpr std::string_view YELLOW  = "\033[33m";
  static constexpr std::string_view BLUE    = "\033[34m";
  static constexpr std::string_view MAGENTA = "\033[35m";
  static constexpr std::string_view CYAN    = "\033[36m";
  static constexpr std::string_view WHITE   = "\033[37m";
  static constexpr std::string_view DEFAULT = "\033[39m";

  static constexpr std::string_view BRIGHT_BLACK   = "\033[90m";
  static constexpr std::string_view BRIGHT_RED     = "\033[91m";
  static constexpr std::string_view BRIGHT_GREEN   = "\033[92m";
  static constexpr std::string_view BRIGHT_YELLOW  = "\033[93m";
  static constexpr std::string_view BRIGHT_BLUE    = "\033[94m";
  static constexpr std::string_view BRIGHT_MAGENTA = "\033[95m";
  static constexpr std::string_view BRIGHT_CYAN    = "\033[96m";
  static constexpr std::string_view BRIGHT_WHITE   = "\033[97m";

  static constexpr std::string_view BG_BLACK   = "\033[40m";
  static constexpr std::string_view BG_RED     = "\033[41m";
  static constexpr std::string_view BG_GREEN   = "\033[42m";
  static constexpr std::string_view BG_YELLOW  = "\033[43m";
  static constexpr std::string_view BG_BLUE    = "\033[44m";
  static constexpr std::string_view BG_MAGENTA = "\033[45m";
  static constexpr std::string_view BG_CYAN    = "\033[46m";
  static constexpr std::string_view BG_WHITE   = "\033[47m";
  static constexpr std::string_view BG_DEFAULT = "\033[49m";

  static constexpr std::string_view BG_BRIGHT_BLACK   = "\033[100m";
  static constexpr std::string_view BG_BRIGHT_RED     = "\033[101m";
  static constexpr std::string_view BG_BRIGHT_GREEN   = "\033[102m";
  static constexpr std::string_view BG_BRIGHT_YELLOW  = "\033[103m";
  static constexpr std::string_view BG_BRIGHT_BLUE    = "\033[104m";
  static constexpr std::string_view BG_BRIGHT_MAGENTA = "\033[105m";
  static constexpr std::string_view BG_BRIGHT_CYAN    = "\033[106m";
  static constexpr std::string_view BG_BRIGHT_WHITE   = "\033[107m";

  static constexpr std::string_view MY_ORANGE = "\033[38;5;166m";
  static constexpr std::string_view MY_PURPLE = "\033[38;5;91m";
  static constexpr std::string_view MY_BLUE   = "\033[38;5;33m";
  static constexpr std::string_view MY_GRAY   = "\033[38;5;240m";

  // NOLINTEND(readability-identifier-naming)

  /**
   * @brief Wrap text in a single style.
   */
  static std::string
  colorize(const StringLike auto &text, std::string_view code)
  {
    return std::string(code) + std::string(text) + std::string(RESET);
  }

  /**
   * @brief Wrap text in multiple styles.
   */
  static std::string
  colorize(const StringLike auto &text, std::initializer_list<std::string_view> codes)
  {
    std::string prefix;
    for (auto code : codes)
      {
        prefix += code;
      }
    return prefix + std::string(text) + std::string(RESET);
  }

  /**
   * @brief Wrap text in multiple styles.
   */
  template <typename... Codes>
  requires(std::convertible_to<Codes, std::string_view> && ...)
  static std::string
  colorize(const StringLike auto &text, [[maybe_unused]] Codes &&...codes)
  {
    std::string prefix;
    ((prefix += std::string_view(codes)), ...);
    return prefix + std::string(text) + std::string(RESET);
  }

  /**
   * @brief Tessellate text with given styles.
   */
  static std::string
  tessellate(const StringLike auto                  &text,
             std::initializer_list<std::string_view> codes,
             std::size_t                             block_size = 1,
             std::string_view                        modifier   = "")
  {
    return tessellate_impl(std::string_view(text),
                           std::span(codes.begin(), codes.size()),
                           block_size,
                           modifier);
  }

  /**
   * @brief Tessellate text with given styles.
   */
  template <typename... Codes>
  requires(std::convertible_to<Codes, std::string_view> && ...)
  static std::string
  tessellate(const StringLike auto &text,
             std::size_t            block_size,
             std::string_view       modifier,
             [[maybe_unused]] Codes &&...codes)
  {
    auto code_arr = std::array {std::string_view(codes)...};
    return tessellate_impl(std::string_view(text),
                           std::span(code_arr),
                           block_size,
                           modifier);
  }

  /**
   * @brief Whether terminal supports ANSI codes.
   */
  [[nodiscard]] static bool
  is_supported() noexcept
  {
    return std::getenv("TERM") != nullptr || std::getenv("COLORTERM") != nullptr;
  }

private:
  static std::string
  tessellate_impl(std::string_view                  text,
                  std::span<const std::string_view> codes,
                  std::size_t                       block_size,
                  std::string_view                  modifier = "")
  {
    if (codes.empty() || text.empty())
      {
        return std::string(text);
      }

    static std::mt19937                        rng {std::random_device {}()};
    std::uniform_int_distribution<std::size_t> dist(0, codes.size() - 1);

    std::string result;

    std::size_t index = 0;
    while (index < text.size())
      {
        result += modifier;
        std::string_view code = codes[dist(rng)];
        result += code;

        // We deviate from the block size a little because we need to use the
        // UTF-8 characters
        std::size_t chars_written = 0;
        std::size_t next_index    = index;
        while (next_index < text.size() && chars_written < block_size)
          {
            if (text[next_index] == '\n')
              {
                break;
              }

            constexpr uint8_t utf8_continuation_mask    = 0xC0;
            constexpr uint8_t utf8_continuation_pattern = 0x80;
            next_index++;
            while (next_index < text.size() &&
                   (static_cast<uint8_t>(text[next_index]) & utf8_continuation_mask) ==
                     utf8_continuation_pattern)
              {
                next_index++;
              }
            chars_written++;
          }

        result.append(text, index, next_index - index);
        result += RESET;

        if (next_index < text.size() && text[next_index] == '\n')
          {
            result += '\n';
            next_index++;
          }

        index = next_index;
      }

    return result;
  }
};

PRISMS_PF_END_NAMESPACE
