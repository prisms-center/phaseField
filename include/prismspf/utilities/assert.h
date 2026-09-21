// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/utilities/logger.h>

#include <prismspf/config.h>

#include <cstdio>
#include <libassert/assert.hpp>
#include <stdexcept>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * There's one reason for this file. The libassert/assert.hpp header may not show up in
 * LSPs without having built the project. This is due to how ExternalProject works and I
 * don't see a good reason around it. LSPs should still be able to autocomplete and
 * include the prismspf/utilities/assert.h header.
 *
 * We include header macro guards for the LSP too
 */
#ifndef DEBUG_ASSERT
#  define DEBUG_ASSERT (void);
#endif

#ifndef ASSERT
#  define ASSERT (void);
#endif

#ifndef ASSUME
#  define ASSUME (void);
#endif

#ifndef PANIC
#  define PANIC (void);
#endif

#ifndef UNREACHABLE
#  define UNREACHABLE(void) ;
#endif

/**
 * We want our own custom failure handler for libassert. There are two things we want to
 * do:
 *   1. Print the assertion to log file
 *   2. Throw an exception for DEBUG_ASSERT and ASSERT rather than abort
 */
[[noreturn]] void
failure_handler(const libassert::assertion_info &info)
{
  libassert::enable_virtual_terminal_processing_if_needed();

  std::string message =
    info.to_string(libassert::terminal_width(libassert::stderr_fileno),
                   libassert::isatty(libassert::stderr_fileno)
                     ? libassert::get_color_scheme()
                     : libassert::color_scheme::blank);
  std::cerr << message << std::endl;

  switch (info.type)
    {
      case libassert::assert_type::assertion:
      case libassert::assert_type::debug_assertion:
        throw std::runtime_error(message);
      case libassert::assert_type::assumption:
      case libassert::assert_type::panic:
      case libassert::assert_type::unreachable:
        std::fflush(stderr);
        std::abort();
      default:
        std::cerr << "Critical error: Unknown libassert::assert_type" << std::endl;
        std::abort();
    }
}

PRISMS_PF_END_NAMESPACE
