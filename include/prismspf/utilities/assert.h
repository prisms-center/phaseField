#pragma once

#include <prismspf/config.h>

#include <libassert/assert.hpp>

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

PRISMS_PF_END_NAMESPACE
