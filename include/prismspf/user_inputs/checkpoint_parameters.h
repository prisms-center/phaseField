// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef checkpoint_parameters_h
#define checkpoint_parameters_h

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/user_inputs/temporal_discretization.h>
#include <prismspf/utilities/utilities.h>

#include <climits>
#include <set>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Struct that holds checkpoint parameters.
 */
struct checkpointParameters
{
public:
  /**
   * \brief Return is the current increment should be checkpointed.
   */
  [[nodiscard]] bool
  should_checkpoint(unsigned int increment) const;

  /**
   * \brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate(const temporalDiscretization &temporal_discretization);

  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Whether to load from a checkpoint
  bool load_from_checkpoint = false;

  // Checkpoint condition type
  std::string condition;

  // Number of checkpoints
  unsigned int n_checkpoints = 0;

  // User given checkpoint list
  std::vector<int> user_checkpoint_list;

  // List of increments for checkpoints
  std::set<unsigned int> checkpoint_list;
};

inline bool
checkpointParameters::should_checkpoint(unsigned int increment) const
{
  return checkpoint_list.find(increment) != checkpoint_list.end();
}

inline void
checkpointParameters::postprocess_and_validate(
  const temporalDiscretization &temporal_discretization)
{
  // If the user has specified a list and we have list checkpoint use that and return
  // early
  if (condition == "LIST")
    {
      for (const auto &increment : user_checkpoint_list)
        {
          checkpoint_list.insert(static_cast<unsigned int>(increment));
        }
      return;
    }

  // If the number of checkpoints is 0 return early
  if (n_checkpoints == 0)
    {
      return;
    }

  // If the number of outputs is greater than the number of increments, force them to be
  // equivalent
  n_checkpoints = std::min(n_checkpoints, temporal_discretization.total_increments);

  // Determine the output list from the other criteria
  if (condition == "EQUAL_SPACING")
    {
      for (unsigned int iteration = 0;
           iteration <= temporal_discretization.total_increments;
           iteration += temporal_discretization.total_increments / n_checkpoints)
        {
          checkpoint_list.insert(iteration);
        }
    }
  else if (condition == "LOG_SPACING")
    {
      checkpoint_list.insert(0);
      for (unsigned int output = 1; output <= n_checkpoints; output++)
        {
          checkpoint_list.insert(static_cast<unsigned int>(std::round(
            std::pow(static_cast<double>(temporal_discretization.total_increments),
                     static_cast<double>(output) / static_cast<double>(n_checkpoints)))));
        }
    }
  else if (condition == "N_PER_DECADE")
    {
      AssertThrow(temporal_discretization.total_increments > 1,
                  dealii::ExcMessage("For n per decaded spaced outputs, the number of "
                                     "increments must be greater than 1."));

      checkpoint_list.insert(0);
      checkpoint_list.insert(1);
      for (unsigned int iteration = 2;
           iteration <= temporal_discretization.total_increments;
           iteration++)
        {
          const auto decade = static_cast<unsigned int>(std::ceil(std::log10(iteration)));
          const auto step_size =
            static_cast<unsigned int>(std::pow(10, decade) / n_checkpoints);
          if (iteration % step_size == 0)
            {
              checkpoint_list.insert(iteration);
            }
        }
    }
  else
    {
      AssertThrow(false, UnreachableCode());
    }
}

inline void
checkpointParameters::print_parameter_summary() const
{
  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Checkpoint Parameters\n"
    << "================================================\n"
    << "Checkpoint condition: " << condition << "\n"
    << "Number of checkpoints: " << n_checkpoints << "\n";

  conditionalOStreams::pout_summary() << "Checkpoint iteration list: ";
  for (const auto &iteration : checkpoint_list)
    {
      conditionalOStreams::pout_summary() << iteration << " ";
    }
  conditionalOStreams::pout_summary() << "\n\n" << std::flush;
}

PRISMS_PF_END_NAMESPACE

#endif
