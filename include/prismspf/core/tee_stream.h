// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/config.h>

#include <ostream>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Combined output streams so we can output to terminal and the summary.log with a
 * single statement.
 */
class TeeStream : public std::ostream
{
public:
  TeeStream(std::ostream &stream1, std::ostream &stream2)
    : std::ostream(&tee_buffer)
    , tee_buffer(stream1, stream2)
  {}

private:
  class TeeBuffer : public std::streambuf
  {
  public:
    TeeBuffer(std::ostream &stream1, std::ostream &stream2)
      : stream1(stream1.rdbuf())
      , stream2(stream2.rdbuf())
    {}

  protected:
    int
    overflow(int character) override
    {
      if (character == EOF)
        {
          return static_cast<int>(!EOF);
        }

      const int result1 = stream1->sputc(static_cast<char_type>(character));
      const int result2 = stream2->sputc(static_cast<char_type>(character));
      return result1 == EOF || result2 == EOF ? EOF : character;
    }

    int
    sync() override
    {
      const int result1 = stream1->pubsync();
      const int result2 = stream2->pubsync();
      return result1 == 0 && result2 == 0 ? 0 : -1;
    }

  private:
    std::streambuf *stream1;
    std::streambuf *stream2;
  };

  TeeBuffer tee_buffer;
};

PRISMS_PF_END_NAMESPACE
