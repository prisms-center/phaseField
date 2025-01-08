#ifndef output_parameters_h
#define output_parameters_h

#include <core/conditional_ostreams.h>
#include <core/user_inputs/temporal_discretization.h>
#include <utilities/utilities.h>

/**
 * \brief Class that holds output parameters.
 */
class outputParameters
{
public:
  /**
   * \brief Constructor.
   */
  outputParameters() = default;

  /**
   * \brief Destructor.
   */
  ~outputParameters() = default;

  /**
   * \brief Compute the list of increments where the solution is output to file.
   */
  void
  compute_output_list(const temporalDiscretization &temporal_discretization);

  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Output file type
  std::string output_file_type;

  // Output file name
  std::string output_file_name;

  // The number of steps between outputting
  uint output_frequency = UINT_MAX;

  // Output condition type
  std::string output_condition;

  // Number of outputs
  uint n_outputs = 0;

  // User given output_list
  std::vector<int> user_output_list;

  // Whether to print timing information with output
  bool print_timing_with_output = false;

  // List of increments that output the solution to file
  std::vector<uint> output_list;
};

inline void
outputParameters::compute_output_list(
  const temporalDiscretization &temporal_discretization)
{
  // If the user has specified a list and we have list output use that and return early
  if (output_condition == "LIST")
    {
      for (const auto &increment : user_output_list)
        {
          output_list.push_back(static_cast<uint>(increment));
        }
      return;
    }

  // If the number of outputs is 0 return early
  if (n_outputs == 0)
    {
      return;
    }

  // If the number of outputs is greater than the number of increments, force them to be
  // equivalent
  n_outputs = std::min(n_outputs, temporal_discretization.total_increments);

  // Determine the output list from the other criteria
  if (output_condition == "EQUAL_SPACING")
    {
      for (uint iteration = 0; iteration <= temporal_discretization.total_increments;
           iteration += temporal_discretization.total_increments / n_outputs)
        {
          output_list.push_back(iteration);
        }
    }
  else if (output_condition == "LOG_SPACING")
    {
      output_list.push_back(0);
      for (uint output = 1; output <= n_outputs; output++)
        {
          output_list.push_back(static_cast<uint>(std::round(
            std::pow(static_cast<double>(temporal_discretization.total_increments),
                     static_cast<double>(output) / static_cast<double>(n_outputs)))));
        }
    }
  else if (output_condition == "N_PER_DECADE")
    {
      AssertThrow(temporal_discretization.total_increments > 1,
                  dealii::ExcMessage("For n per decaded spaced outputs, the number of "
                                     "increments must be greater than 1."));

      output_list.push_back(0);
      output_list.push_back(1);
      for (uint iteration = 2; iteration <= temporal_discretization.total_increments;
           iteration++)
        {
          const uint decade    = std::ceil(std::log10(iteration));
          const auto step_size = static_cast<uint>(std::pow(10, decade) / n_outputs);
          if (iteration % step_size == 0)
            {
              output_list.push_back(iteration);
            }
        }
    }
  else
    {
      AssertThrow(false, dealii::ExcMessage("Invalid output spacing type."));
    }
}

inline void
outputParameters::print_parameter_summary() const
{
  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "\tOutput Parameters\n"
    << "================================================\n"
    << "Output file type: " << output_file_type << "\n"
    << "Output file name: " << output_file_name << "\n"
    << "Output frequency: " << output_frequency << "\n"
    << "Output condition: " << output_condition << "\n"
    << "Number of outputs: " << n_outputs << "\n"
    << "Print timing info: " << bool_to_string(print_timing_with_output) << "\n";

  conditionalOStreams::pout_summary() << "Output iteration list: ";
  for (const auto &iteration : output_list)
    {
      conditionalOStreams::pout_summary() << iteration << " ";
    }
  conditionalOStreams::pout_summary() << "\n\n" << std::flush;
}

#endif