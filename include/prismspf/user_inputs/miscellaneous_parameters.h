// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mpi.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <mpi.h>
#include <random>

PRISMS_PF_BEGIN_NAMESPACE
using RNGEngine = std::mt19937;

/**
 * @brief Struct that holds miscellaneous parameters.
 */
struct MiscellaneousParameters
{
public:
  /**
   * @brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate();

  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  /**
   * @brief Set the random seed and initialize the RNG.
   */
  void
  set_random_seed(const unsigned int &_random_seed)
  {
    random_seed = _random_seed;
    rng.seed(random_seed);
  }

  unsigned int random_seed = 2025;
  // Use a different seed for each MPI process to avoid correlated events
  mutable RNGEngine rng {random_seed};
};

inline void
MiscellaneousParameters::postprocess_and_validate()
{}

inline void
MiscellaneousParameters::print_parameter_summary() const
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Miscellaneous Parameters\n"
    << "================================================\n"
    << "Random seed: " << random_seed << "\n"
    << std::flush;
}

PRISMS_PF_END_NAMESPACE
