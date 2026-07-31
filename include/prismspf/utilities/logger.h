#pragma once

#include <deal.II/base/mpi.h>

#include <prismspf/utilities/assert.h>

#include <prismspf/config.h>

#include <algorithm>
#include <ios>
#include <iostream>
#include <mpi.h>
#include <optional>
#include <ostream>
#include <streambuf>
#include <string>
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

  bool
  is_active() const
  {
    return active;
  }

protected:
  int_type
  overflow(int_type ch) override
  {
    if (traits_type::eq_int_type(ch, traits_type::eof()))
      return traits_type::not_eof(ch);

    if (!active)
      return ch;

    bool ok = true;

    const char c = traits_type::to_char_type(ch);

    for (auto *buffer : buffers)
      {
        if (traits_type::eq_int_type(buffer->sputc(c), traits_type::eof()))
          ok = false;
      }

    return ok ? ch : traits_type::eof();
  }

  std::streamsize
  xsputn(const char *s, std::streamsize n) override
  {
    if (!active)
      return n;

    std::streamsize min_written = n;

    for (auto *buffer : buffers)
      {
        const std::streamsize written = buffer->sputn(s, n);
        min_written                   = std::min(min_written, written);
      }

    return min_written;
  }

  int
  sync() override
  {
    if (!active)
      return 0;

    bool ok = true;

    for (auto *buffer : buffers)
      {
        if (buffer->pubsync() != 0)
          ok = false;
      }

    return ok ? 0 : -1;
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
  explicit LogStream(std::string file)
    : LogStream(std::optional<std::string> {std::move(file)}, std::nullopt)
  {}

  /**
   * Output to std::cout and one file on a given MPI process.
   */
  LogStream(std::string file, unsigned int process_id)
    : LogStream(std::optional<std::string> {std::move(file)},
                std::optional<unsigned int> {process_id})
  {}

  void
  set_condition(const bool condition)
  {
    user_condition = condition;
    fanout.set_active(rank_is_enabled && user_condition);
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
  operator<<(T &&t)
  {
    fanout.get() << std::forward<T>(t);
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
      return true;

    return this_rank() == *process_id;
  }

  static std::optional<std::ofstream>
  open_file_if_needed(const std::optional<std::string> &file,
                      std::optional<unsigned int>       process_id)
  {
    if (!file)
      return std::nullopt;

    if (!rank_matches(process_id))
      return std::nullopt;

    std::string filename = *file;

    // If no specific process was requested, create one file per rank.
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

  LogStream(std::optional<std::string> file, std::optional<unsigned int> process_id)
    : rank_is_enabled(rank_matches(process_id))
    , user_condition(true)
    , file_stream(open_file_if_needed(file, process_id))
  {
    fanout.add_stream(std::cout);

    if (file_stream)
      fanout.add_stream(*file_stream);

    fanout.set_active(rank_is_enabled && user_condition);
  }

private:
  bool                         rank_is_enabled = true;
  bool                         user_condition  = true;
  std::optional<std::ofstream> file_stream;
  FanoutOStream                fanout;
};

/**
 * @brief Format string according to some logging style
 */
class LogFormatter
{};

/**
 * @brief Main logging class
 *
 * This class manages instances of `LogStream` and passes strings through the formatting
 * utilities to keep a consistent style.
 */
class Logger
{};

PRISMS_PF_END_NAMESPACE
