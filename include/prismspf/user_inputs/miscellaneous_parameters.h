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
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override;
  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) override;

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override;
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
  set_random_seed(const unsigned int &_random_seed);

  // Random seed
  unsigned int random_seed = 2025;

  // RNG
  mutable RNGEngine rng {random_seed};
};

PRISMS_PF_END_NAMESPACE
