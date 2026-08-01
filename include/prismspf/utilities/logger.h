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
  add_stream(std::ostream &out)
  {
    buffers.push_back(out.rdbuf());
  }

  void
  set_active(const bool active)
  {
    this->active = active;
  }

  [[nodiscard]] bool
  is_active() const
  {
    return active;
  }

protected:
  int_type
  overflow(int_type int_char) override
  {
    if (traits_type::eq_int_type(int_char, traits_type::eof()))
      {
        return traits_type::not_eof(int_char);
      }
    if (!active)
      {
        return int_char;
      }
    bool       is_ok = true;
    const char _char = traits_type::to_char_type(int_char);
    for (auto *buffer : buffers)
      {
        if (traits_type::eq_int_type(buffer->sputc(_char), traits_type::eof()))
          {
            is_ok = false;
          }
      }
    return is_ok ? int_char : traits_type::eof();
  }

  std::streamsize
  xsputn(const char *_char, std::streamsize n) override
  {
    if (!active)
      {
        return n;
      }
    std::streamsize min_written = n;
    for (auto *buffer : buffers)
      {
        const std::streamsize written = buffer->sputn(_char, n);
        min_written                   = std::min(min_written, written);
      }
    return min_written;
  }

  int
  sync() override
  {
    if (!active)
      {
        return 0;
      }
    bool is_ok = true;
    for (auto *buffer : buffers)
      {
        if (buffer->pubsync() != 0)
          {
            is_ok = false;
          }
      }
    return is_ok ? 0 : -1;
  }

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
  FanoutOStream()
    : stream(&buffer)
  {}

  void
  add_stream(std::ostream &out)
  {
    buffer.add_stream(out);
  }

  void
  set_active(const bool active)
  {
    buffer.set_active(active);
  }

  bool
  is_active() const
  {
    return buffer.is_active();
  }

  std::ostream &
  get()
  {
    return stream;
  }

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
  /**
   * Output to std::cout on all MPI processes.
   */
  LogStream()
    : LogStream(std::nullopt, std::nullopt)
  {}

  /**
   * Output to std::cout on a given MPI process.
   */
  explicit LogStream(unsigned int process_id)
    : LogStream(std::nullopt, std::optional<unsigned int> {process_id})
  {}

  /**
   * Output to std::cout and one file per MPI process.
   *
   * The file name is appended by the MPI process id.
   */
  explicit LogStream(const std::string &file)
    : LogStream(std::optional<std::string> {file}, std::nullopt)
  {}

  /**
   * Output to std::cout and one file on a given MPI process.
   */
  LogStream(const std::string &file, unsigned int process_id)
    : LogStream(std::optional<std::string> {file},
                std::optional<unsigned int> {process_id})
  {}

  void
  add_file([[maybe_unused]] const std::string &file,
           [[maybe_unused]] unsigned int       process_id)
  {
    if (file_stream)
      {
        return;
      }
    file_stream = open_file_if_needed(file, process_id);
  }

  void
  set_condition(const bool condition)
  {
    user_condition = condition;
    fanout.set_active(is_active());
  }

  bool
  is_active() const
  {
    return rank_is_enabled && user_condition;
  }

  std::ostream &
  get()
  {
    return fanout.get();
  }

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
  static unsigned int
  this_rank()
  {
    return dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  }

  static bool
  rank_matches(std::optional<unsigned int> process_id)
  {
    if (!process_id)
      {
        return true;
      }
    return this_rank() == *process_id;
  }

  static std::optional<std::ofstream>
  open_file_if_needed(const std::optional<std::string> &file,
                      std::optional<unsigned int>       process_id)
  {
    if (!file)
      {
        return std::nullopt;
      }
    if (!rank_matches(process_id))
      {
        return std::nullopt;
      }

    // If no specific process was requested, create one file per rank.
    std::string filename = *file;
    if (!process_id)
      {
        filename += "_" + std::to_string(this_rank());
      }
    std::optional<std::ofstream> out {std::in_place,
                                      filename,
                                      std::ios::out | std::ios::trunc};
    ASSERT(*out, "Could not open log file", filename);
    return out;
  }

  LogStream(const std::optional<std::string> &file,
            std::optional<unsigned int>       process_id)
    : rank_is_enabled(rank_matches(process_id))
    , file_stream(open_file_if_needed(file, process_id))
  {
    fanout.add_stream(std::cout);
    if (file_stream)
      {
        fanout.add_stream(*file_stream);
      }
    fanout.set_active(rank_is_enabled && user_condition);
  }

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
    RainbowTessellate // Rainbow tesselation
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

  friend std::ostream &
  operator<<(std::ostream &stream, const Manipulator &manip)
  {
    return stream << format(manip.section, manip.style, manip.text);
  }

private:
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
    switch (_section)
      {
        case Section::Normal:
          return std::string {string};
        case Section::Title:
          return std::string {PRISMS_PF_ASCII};
        case Section::Subtitle:
          return "version " + std::string {PRISMS_PF_VERSION} + "\n" + "git revision " +
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
  format(Section _section, Style _style, const std::string_view &text);
};

/**
 * @brief Main logging class
 *
 * This class manages two static instances of the LogStream class: one for std::cout and
 * one fanout stream for std::cout and some log file. Both of these only run on the 0th
 * process.
 *
 * We must also consider that certain information is only relevant when certain CMake
 * flags are on (e.g., DEBUG and VERBOSE). This produces a lot of separate ostreams that
 * would have to be managed in the code. Rather than do that, we rely on manipulators to
 * tell us what strings go where and how they should be formatted.
 *
 * Thus, this code here:
 *
 * @code{.cpp}
 * @endcode
 *
 * Would produce the following resuts:
 *
 */
class Logger
{
public:
  /**
   * @brief Clear terminal and print logo upon construction
   */
  Logger()
  {
    // NOTE: This won't print to the log file because the class gets constructed before it
    // is set. Since it's just the logo it doesn't matter.
    if (colorize())
      {
        conditional_file_stream << TerminalColor::ERASE_SCREEN;
      }
    conditional_file_stream << LogFormatter::title() << LogFormatter::subtitle()
                            << std::flush;
  }

  /**
   * @brief Set the log file to write to.
   */
  static void
  set_file(const std::string &file)
  {
    conditional_file_stream.add_file(file, 0);
  }

  /**
   * @brief Static instance
   */
  static Logger &
  instance()
  {
    static Logger logger;
    return logger;
  }

  /**
   * @brief Stream operator
   */
  template <typename T>
  Logger &
  operator<<(const T &type)
  {
    conditional_stream << type;
    return *this;
  }

  /**
   * @brief Stream operator
   */
  Logger &
  operator<<(std::ostream &(*manip)(std::ostream &) )
  {
    conditional_stream << manip;
    return *this;
  }

  /**
   * @brief Whether the terminal supports color
   */
  [[nodiscard]] static bool
  colorize()
  {
    return term_supports_color;
  }

private:
  inline static bool         term_supports_color {TerminalColor::is_supported()};
  inline static unsigned int indent_level = 0;
  inline static LogStream    conditional_stream {0};
  inline static LogStream    conditional_file_stream {0};
};

// TODO: Move this
inline std::string
LogFormatter::format(Section _section, Style _style, const std::string_view &text)
{
  if (Logger::colorize())
    {
      return style(section(text, _section), _style);
    }
  return section(text, _section);
}

PRISMS_PF_END_NAMESPACE
