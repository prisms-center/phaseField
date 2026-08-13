#include <prismspf/utilities/logger.h>

#include <prismspf/config.h>

#include "prismspf/utilities/terminal.h"

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

void
FanoutStreamBuf::add_stream(std::ostream &out)
{
  buffers.push_back(out.rdbuf());
}

void
FanoutStreamBuf::set_active(bool _active)
{
  active = _active;
}

bool
FanoutStreamBuf::is_active() const
{
  return active;
}

typename FanoutStreamBuf::int_type
FanoutStreamBuf::overflow(int_type int_char)
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
FanoutStreamBuf::xsputn(const char *_char, std::streamsize n)
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
FanoutStreamBuf::sync()
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

FanoutOStream::FanoutOStream()
  : stream(&buffer)
{}

void
FanoutOStream::add_stream(std::ostream &out)
{
  buffer.add_stream(out);
}

void
FanoutOStream::set_active(bool active)
{
  buffer.set_active(active);
}

bool
FanoutOStream::is_active() const
{
  return buffer.is_active();
}

std::ostream &
FanoutOStream::get()
{
  return stream;
}

LogStream::LogStream(BlankOutputTag tag)
  : LogStream(std::nullopt, std::nullopt, false)
{}

LogStream::LogStream()
  : LogStream(std::nullopt, std::nullopt, true)
{}

LogStream::LogStream(unsigned int process_id)
  : LogStream(std::nullopt, std::optional<unsigned int> {process_id}, true)
{}

LogStream::LogStream(const std::string &file)
  : LogStream(std::optional<std::string> {file}, std::nullopt, true)
{}

LogStream::LogStream(const std::string &file, unsigned int process_id)
  : LogStream(std::optional<std::string> {file},
              std::optional<unsigned int> {process_id},
              true)
{}

void
LogStream::add_file(const std::string &file, unsigned int process_id)
{
  if (file_stream)
    {
      return;
    }
  file_stream = open_file_if_needed(file, process_id);
  if (file_stream)
    {
      fanout.add_stream(*file_stream);
    }
}

void
LogStream::set_condition(bool condition)
{
  user_condition = condition;
  fanout.set_active(is_active());
}

bool
LogStream::is_active() const
{
  return rank_is_enabled && user_condition;
}

std::ostream &
LogStream::get()
{
  return fanout.get();
}

unsigned int
LogStream::this_rank() const
{
  return dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
}

bool
LogStream::rank_matches(std::optional<unsigned int> process_id) const
{
  if (!process_id)
    {
      return true;
    }
  return this_rank() == *process_id;
}

std::optional<std::ofstream>
LogStream::open_file_if_needed(const std::optional<std::string> &file,
                               std::optional<unsigned int>       process_id) const
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

LogStream::LogStream(const std::optional<std::string> &file,
                     std::optional<unsigned int>       process_id,
                     bool                              enable_cout)
  : rank_is_enabled(rank_matches(process_id))
  , file_stream(open_file_if_needed(file, process_id))
{
  if (enable_cout)
    {
      fanout.add_stream(std::cout);
    }
  if (file_stream)
    {
      fanout.add_stream(*file_stream);
    }
  fanout.set_active(rank_is_enabled && user_condition);
}

std::string
LogFormatter::format(Section                 _section,
                     Style                   _style,
                     const std::string_view &text,
                     bool                    use_ansi)
{
  const static bool term_supports_color = TerminalColor::is_supported();
  if (term_supports_color && use_ansi)
    {
      return style(section(text, _section), _style);
    }
  return section(text, _section);
}

Logger::IndentScope::IndentScope()
{
  increment_indent();
}

Logger::IndentScope::~IndentScope()
{
  decrement_indent();
}

Logger::Logger()
{
  // NOTE: We must check that the terminal supports color codes before printing with
  // them to cout
  cout << LogFormatter::title() << LogFormatter::subtitle() << std::flush;
}

void
Logger::set_file(const std::string &file)
{
  instance().log_file.add_file(file, 0);
}

Logger &
Logger::instance()
{
  static Logger logger;
  return logger;
}

void
Logger::increment_indent()
{
  ++instance().indent_level;
}

void
Logger::decrement_indent()
{
  if (instance().indent_level > 0)
    {
      --instance().indent_level;
    }
}

void
Logger::reset_indent()
{
  instance().indent_level = 0;
}

void
Logger::write_indent()
{
  instance().cout << std::string(2 * instance().indent_level, ' ');
  instance().log_file << std::string(2 * instance().indent_level, ' ');
}

PRISMS_PF_END_NAMESPACE
