// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>

#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <exception>

PRISMS_PF_BEGIN_NAMESPACE

// NOLINTBEGIN(readability-identifier-naming)

/**
 * Macro for assertion that is only called in DEBUG mode. This is used to
 * throw a catchable exception with zero overhead in RELEASE mode.
 * dealii Assert will call abort(), which is not useful to us.
 */
#ifdef DEBUG
#  define AssertThrowDebug(cond, exc) AssertThrow(cond, exc);
#else
#  define AssertThrowDebug(cond, exc)
#endif

// NOLINTEND(readability-identifier-naming)

// NOLINTBEGIN

class ExcNaN : public std::exception
{
private:
  std::string message;

public:
  // Constructor to store a custom message
  explicit ExcNaN(const std::string &msg = "")
    : message("PRISMS-PF detected somewhere that a value was NaN.\n" + msg)
  {}

  // Overriding the virtual what() method
  const char *
  what() const noexcept override
  {
    return message.c_str();
  }
};

/**
 * Exception for parts of the library that have yet to be implemented yet. The argument is
 * used to provide additional context for the feature that has yet to be implemented.
 */
DeclException1(
  FeatureNotImplemented,
  std::string,
  << "The following feature has yet to be implemented in PRISMS-PF:\n  " << arg1
  << "\nCheck the issues section of PRISMS-PF's github to see if this feature is under "
     "development. Additionally, please considering provided a patch to PRISMS-PF if you "
     "feel that feature is worthwhile for yourself and others.");

/**
 * Exception for when we have reached part of the code that should not have been reached.
 * This is common in switch statements and conditional chains.
 */
DeclExceptionMsg(UnreachableCode, "This code should not have been reached.");

/**
 * Exception for a user trying to access a variable in from VariableContainer that has not
 * been specified as a dependency.
 */
DeclException2(DependencyNotFound,
               Types::Index,
               std::string,
               << "Attempted access of the variable with index " << arg1
               << " and dependency type " << arg2
               << " that was not marked as needed. Please check CustomAttributeLoader.");

// NOLINTEND

PRISMS_PF_END_NAMESPACE
