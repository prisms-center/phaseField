#pragma once

#include <deal.II/base/mpi.h>

#include <prismspf/utilities/assert.h>
#include <prismspf/utilities/terminal.h>

#include <prismspf/config.h>

#include <algorithm>
#include <fstream>
#include <ios>
#include <iostream>
#include <mpi.h>
#include <optional>
#include <ostream>
#include <streambuf>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * Fanout stream buffer
 */
class FanoutStreamBuf : public std::streambuf
{
public:
  void
  add_stream(std::ostream &out);

  void
  set_active(bool _active);

  [[nodiscard]] bool
  is_active() const;

protected:
  int_type
  overflow(int_type int_char) override;

  std::streamsize
  xsputn(const char *_char, std::streamsize n) override;

  int
  sync() override;

private:
  bool                          active = true;
  std::vector<std::streambuf *> buffers;
};

/**
 * Fanout output stream
 */
class FanoutOStream
{
public:
  FanoutOStream();

  void
  add_stream(std::ostream &out);

  void
  set_active(bool active);

  bool
  is_active() const;

  std::ostream &
  get();

private:
  FanoutStreamBuf buffer;
  std::ostream    stream;
};

/**
 * @brief Logging of data to a given ostream
 */
class LogStream
{
public:
  struct BlankOutputTag
  {};

  /**
   * Blank output
   */
  explicit LogStream([[maybe_unused]] BlankOutputTag tag);

  /**
   * Output to std::cout on all MPI processes.
   */
  LogStream();

  /**
   * Output to std::cout on a given MPI process.
   */
  explicit LogStream(unsigned int process_id);

  /**
   * Output to std::cout and one file per MPI process.
   *
   * The file name is appended by the MPI process id.
   */
  explicit LogStream(const std::string &file);

  /**
   * Output to std::cout and one file on a given MPI process.
   */
  LogStream(const std::string &file, unsigned int process_id);

  void
  add_file([[maybe_unused]] const std::string &file,
           [[maybe_unused]] unsigned int       process_id);

  void
  set_condition(bool condition);

  bool
  is_active() const;

  std::ostream &
  get();

  template <typename T>
  LogStream &
  operator<<(T &&type)
  {
    fanout.get() << std::forward<T>(type);
    return *this;
  }

  LogStream &
  operator<<(std::ostream &(*manip)(std::ostream &) )
  {
    fanout.get() << manip;
    return *this;
  }

private:
  unsigned int
  this_rank() const;

  bool
  rank_matches(std::optional<unsigned int> process_id) const;

  std::optional<std::ofstream>
  open_file_if_needed(const std::optional<std::string> &file,
                      std::optional<unsigned int>       process_id) const;

  LogStream(const std::optional<std::string> &file,
            std::optional<unsigned int>       process_id,
            bool                              enable_cout);

  bool                         rank_is_enabled = true;
  bool                         user_condition  = true;
  std::optional<std::ofstream> file_stream;
  FanoutOStream                fanout;
};

/**
 * @brief Format string according to some logging style
 *
 * Primarily relies on the TerminalColor colorize functions.
 */
class LogFormatter
{
public:
  /**
   * Constructor
   *
   * Deleted because the class is only intended to be used with stream operators.
   */
  LogFormatter() = delete;

  /**
   * Available styles.
   */
  enum class Style
  {
    Normal,           // None
    Info,             // Blue
    Warning,          // Red
    Error,            // Yellow
    Success,          // Green
    Debug,            // Grey
    Verbose,          // Grey
    RainbowTessellate // Rainbow tessellation
  };

  /**
   * Available section styles.
   */
  enum class Section
  {
    Normal,    // None
    Title,     // Tessellated ASCII art banner
    Subtitle,  // Version info
    Section,   // ===== Section Name =====
    Subsection // --- Subsection Name ---
  };

  /**
   * Manipulator struct
   */
  struct Manipulator
  {
    Section     section = Section::Normal;
    Style       style   = Style::Normal;
    std::string text;
  };

  /**
   * @brief Format text as a title heading.
   */
  static Manipulator
  title()
  {
    return {Section::Title, Style::RainbowTessellate, ""};
  }

  /**
   * @brief Format text as a subtitle heading.
   */
  static Manipulator
  subtitle()
  {
    return {Section::Subtitle, Style::Info, ""};
  }

  /**
   * @brief Format text as a section heading.
   */
  static Manipulator
  section(std::string string)
  {
    return {Section::Section, Style::Normal, std::move(string)};
  }

  /**
   * @brief Format text as a subsection heading.
   */
  static Manipulator
  subsection(std::string string)
  {
    return {Section::Subsection, Style::Normal, std::move(string)};
  }

  /**
   * @brief Format text as info.
   */
  static Manipulator
  info(std::string string)
  {
    return {Section::Normal, Style::Info, std::move(string)};
  }

  /**
   * @brief Format text as warning.
   */
  static Manipulator
  warning(std::string string)
  {
    return {Section::Normal, Style::Warning, std::move(string)};
  }

  /**
   * @brief Format text as error.
   */
  static Manipulator
  error(std::string string)
  {
    return {Section::Normal, Style::Error, std::move(string)};
  }

  /**
   * @brief Format text as success.
   */
  static Manipulator
  success(std::string string)
  {
    return {Section::Normal, Style::Success, std::move(string)};
  }

  /**
   * @brief Format text as debug.
   */
  static Manipulator
  debug(std::string string)
  {
    return {Section::Normal, Style::Debug, std::move(string)};
  }

  static std::string
  style(const std::string_view &string, Style _style)
  {
    switch (_style)
      {
        case Style::Normal:
          return TerminalColor::colorize(string, TerminalColor::RESET);
        case Style::Info:
          return TerminalColor::colorize(string, TerminalColor::CYAN);
        case Style::Warning:
          return TerminalColor::colorize(string, TerminalColor::YELLOW);
        case Style::Error:
          return TerminalColor::colorize(string,
                                         {TerminalColor::BOLD,
                                          TerminalColor::BRIGHT_RED});
        case Style::Success:
          return TerminalColor::colorize(string, TerminalColor::GREEN);
        case Style::Debug:
        case Style::Verbose:
          return TerminalColor::colorize(string,
                                         {TerminalColor::DIM,
                                          TerminalColor::BRIGHT_BLACK});
        case Style::RainbowTessellate:
          return TerminalColor::tessellate(string,
                                           {TerminalColor::BRIGHT_RED,
                                            TerminalColor::BRIGHT_GREEN,
                                            TerminalColor::BRIGHT_YELLOW,
                                            TerminalColor::BRIGHT_BLUE,
                                            TerminalColor::BRIGHT_MAGENTA,
                                            TerminalColor::BRIGHT_CYAN,
                                            TerminalColor::BRIGHT_WHITE},
                                           2);
        default:
          return TerminalColor::colorize(string, TerminalColor::RESET);
      }
  }

  static std::string
  section(const std::string_view &string, Section _section)
  {
    // TODO: Avoid the 79 magic number here
    switch (_section)
      {
        case Section::Normal:
          return std::string {string};
        case Section::Title:
          return std::string {PRISMS_PF_ASCII};
        case Section::Subtitle:
          return "Version " + std::string {PRISMS_PF_VERSION} + "\n" + "Git revision " +
                 std::string {PRISMS_PF_GIT_VERSION} + "\n";
        case Section::Section:
          return "\n" + std::string(79, '=') + "\n" + std::string {string} + "\n" +
                 std::string(79, '=') + "\n";
        case Section::Subsection:
          return "\n" + std::string(79, '-') + "\n" + std::string {string} + "\n" +
                 std::string(79, '-') + "\n";
        default:
          return std::string {string};
      }
  }

  static std::string
  format(Section                 _section,
         Style                   _style,
         const std::string_view &text,
         bool                    use_ansi = true);

  friend std::ostream &
  operator<<(std::ostream &stream, const LogFormatter::Manipulator &manip)
  {
    return stream << LogFormatter::format(manip.section, manip.style, manip.text);
  }
};

/**
 * @brief Main logging class
 *
 * This class manages two static instances of the LogStream class: one for std::cout one
 * for some log file. Both of these only run on the 0th process.
 *
 * We must also consider that certain information is only relevant when certain CMake
 * flags are on (e.g., DEBUG and VERBOSE). This produces a lot of separate ostreams that
 * would have to be managed in the code. Rather than do that, we rely on manipulators to
 * tell us what strings go where and how they should be formatted.
 *
 * @note Indents are only applied when starting on a newline. For now, this is only
 * detected after using `std::endl`.
 *
 * @warn When using the verbose style do not follow with `std::endl`. This is handled
 * uniquely by that style so as to not produce extra whitespace.
 *
 * Thus, this code here:
 *
 * @code{.cpp}
 * @endcode
 *
 * Would produce the following results:
 *
 */
class Logger : public std::ostream
{
public:
  /**
   * @brief Scope indents
   */
  class IndentScope
  {
  public:
    IndentScope();
    ~IndentScope();
  };

  /**
   * @brief Clear terminal and print logo upon construction
   */
  Logger();

  /**
   * @brief Set the log file to write to.
   */
  static void
  set_file(const std::string &file);

  /**
   * @brief Static instance
   *
   * Use this to access the logger in the code.
   */
  static Logger &
  instance();

  /**
   * @brief Increment indentation level
   */
  static void
  increment_indent();

  /**
   * @brief Decrement indentation level
   */
  static void
  decrement_indent();

  /**
   * @brief Reset indentation level
   *
   * @note This gets called when you use a section or subsection manipulator
   */
  static void
  reset_indent();

  /**
   * @brief Stream operator
   */
  Logger &
  operator<<(const LogFormatter::Manipulator &manip)
  {
    if (manip.section == LogFormatter::Section::Section ||
        manip.section == LogFormatter::Section::Subsection)
      {
        reset_indent();
      }

    if (at_line_start)
      {
        write_indent();
        at_line_start = false;
      }

    log_file << LogFormatter::format(manip.section, manip.style, manip.text, false);
    if (manip.style == LogFormatter::Style::Verbose)
      {
        log_file << std::endl;
        at_line_start = true;
      }
    else
      {
        cout << LogFormatter::format(manip.section, manip.style, manip.text);
      }

    return *this;
  }

  /**
   * @brief Stream operator
   */
  template <typename T>
  Logger &
  operator<<(const T &type)
  {
    if (at_line_start)
      {
        write_indent();
        at_line_start = false;
      }

    cout << type;
    log_file << type;
    return *this;
  }

  /**
   * @brief Stream operator
   */
  Logger &
  operator<<(std::ostream &(*manip)(std::ostream &) )
  {
    cout << manip;
    log_file << manip;
    at_line_start = true;
    return *this;
  }

private:
  static void
  write_indent();

  unsigned int indent_level  = 0;
  bool         at_line_start = true;

  LogStream cout {0};
  LogStream log_file {LogStream::BlankOutputTag {}};
};

PRISMS_PF_END_NAMESPACE
