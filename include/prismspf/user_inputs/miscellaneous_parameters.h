// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>

#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <mpi.h>
#include <random>

PRISMS_PF_BEGIN_NAMESPACE

using RNGEngine = std::mt19937;

/**
 * @brief Struct that holds miscellaneous parameters.
 */
struct MiscellaneousParameters : public ParameterBase
{
public:
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override
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
  };

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) override
  {
    parameter_handler.enter_subsection("miscellaneous");
    {
      set_random_seed((unsigned int) (parameter_handler.get_integer("random seed")));
    }
    parameter_handler.leave_subsection();
  };

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override {
    // TODO: Do this later
  };

  /**
   * @brief Set the random seed and initialize the RNG.
   *
   * NOTE: When using this method to set the RNG, the MPI process number is added to the
   * random seed to avoid correlated events between processors. This means that to
   * reproduce the same results as someone else you need both the number of MPI processes
   * and the random seed.
   *
   * If you would like to change the default behavior you can override the rng manually.
   */
  void
  set_random_seed(const unsigned int &_random_seed)
  {
    random_seed = _random_seed;
    rng.seed(random_seed + dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
  }

  // Random seed
  unsigned int random_seed = 2025;

  // RNG
  mutable RNGEngine rng {random_seed};
};

PRISMS_PF_END_NAMESPACE
