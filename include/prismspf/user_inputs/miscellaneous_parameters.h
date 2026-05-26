// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

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

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler) const;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler);

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

inline void
MiscellaneousParameters::declare_parameters(
  dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.enter_subsection("miscellaneous");
  {
    parameter_handler.declare_entry(
      "random seed",
      "2025",
      dealii::Patterns::Integer(0, INT_MAX),
      "The random seed for the simulation. "
      "This is used to initialize the random number generator.");
  }
  parameter_handler.leave_subsection();
}

inline void
MiscellaneousParameters::assign_parameters(dealii::ParameterHandler &parameter_handler)
{
  parameter_handler.enter_subsection("miscellaneous");
  {
    set_random_seed(static_cast<unsigned int>(
      parameter_handler.get_integer("random seed") +
      dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)));
  }
  parameter_handler.leave_subsection();
}

PRISMS_PF_END_NAMESPACE
