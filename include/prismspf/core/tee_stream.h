// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef tee_stream_h
#define tee_stream_h

#include <prismspf/config.h>

#include <ostream>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Combined output streams so we can output to terminal and the summary.log with a
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
    overflow(int c) override
    {
      if (c == EOF)
        {
          return static_cast<int>(!EOF);
        }
      else
        {
          const int r1 = stream1->sputc(c);
          const int r2 = stream2->sputc(c);
          return r1 == EOF || r2 == EOF ? EOF : c;
        }
    }

    int
    sync() override
    {
      const int r1 = stream1->pubsync();
      const int r2 = stream2->pubsync();
      return r1 == 0 && r2 == 0 ? 0 : -1;
    }

  private:
    std::streambuf *stream1;
    std::streambuf *stream2;
  };

  TeeBuffer tee_buffer;
};

PRISMS_PF_END_NAMESPACE

#endif