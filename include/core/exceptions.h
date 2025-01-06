#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <deal.II/base/exceptions.h>

using namespace dealii;

/**
 * Exception for parts of the library that have yet to be implemented yet. The argument is
 * used to provide context for the feature that has yet to be implemented.
 */
DeclException1(
  FeatureNotImplemented,
  std::string,
  << "The following feature has yet to be implemented in PRISMS-PF:\n"
  << arg1
  << "\nCheck the issues section of PRISMS-PF's github to see if this feature is under "
     "development. Additionally, please considering provided a patch to PRISMS-PF if you "
     "feel that feature is worthwhile for yourself and others.");

#endif