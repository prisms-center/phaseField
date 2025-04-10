// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mpi.h>

#include <prismspf/core/conditional_ostreams.h>

#include <prismspf/config.h>

#include <fstream>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Class to parse command line options.
 */
class parseCMDOptions
{
public:
  // NOLINTBEGIN (cppcoreguidelines-pro-bounds-pointer-arithmetic)

  parseCMDOptions(int &_argc, char **argv)
    : argc(_argc)
  {
    for (int i = 1; i < argc; ++i)
      {
        tokens.emplace_back(argv[i]);
      }
  }

  // NOLINTEND (cppcoreguidelines-pro-bounds-pointer-arithmetic)

  std::string
  get_parameters_filename()
  {
    // Allowed number of arguments
#ifdef PRISMS_PF_WITH_CALIPER
    const int n_args = 5;
#else
    const int n_args = 3;
#endif

    // Check that there aren't too many arguments
    if (argc > n_args)
      {
        throw std::runtime_error(
          "Too many command line arguments. The arguments should specify "
          "the input file name and caliper configurations, if applicable. For example, "
          "`./main -i parameters.prm`");
      }

    // Default filename
    std::string parameters_filename = "parameters.prm";

    if (argc >= n_args - 1 && !cmd_option_exists("-i"))
      {
        throw std::runtime_error(
          "Invalid command line option. Use `-i` to specify the input file, "
          "e.g., `./main -i parameters.prm`.");
      }

    // If -i is provided, check for an argument
    if (cmd_option_exists("-i"))
      {
        parameters_filename = get_cmd_option("-i");
      }

    // Validate the parameters file (i.e., it ends in .prm)
    if (parameters_filename.size() < 4 ||
        parameters_filename.substr(parameters_filename.size() - 4) != ".prm")
      {
        throw std::runtime_error("The input file must have the `.prm` extension. Please "
                                 "rename the file accordingly.");
      }

    // Check if the file exists
    const std::ifstream ifs_prm(parameters_filename);
    if (!ifs_prm)
      {
        throw std::runtime_error("The specified parameters file `" + parameters_filename +
                                 "` does not exist.");
      }

    // Log the filename being used
    conditionalOStreams::pout_base()
      << "Using the input parameter file: " << parameters_filename << "\n";

    return parameters_filename;
  }

#ifdef PRISMS_PF_WITH_CALIPER
  std::string
  get_caliper_configuration()
  {
    // Allowed number of arguments
    int n_args = 5;

    // Check that there aren't too many arguments
    if (argc > n_args)
      {
        throw std::runtime_error(
          "Too many command line arguments. The arguments should specify "
          "the input file name and caliper configurations, if applicable. For example, "
          "`./main -i parameters.prm -P runtime-report`");
      }

    // Default config
    std::string config = "runtime-report";

    if (argc >= n_args - 1 && !cmd_option_exists("-P"))
      {
        throw std::runtime_error(
          "Invalid command line option. Use `-P` to specify caliper configurations, "
          "e.g., `./main -P runtime-report`.");
      }

    // If -P is provided, check for an argument
    if (cmd_option_exists("-P"))
      {
        config = get_cmd_option("-P");
      }

    // Log the filename being used
    conditionalOStreams::pout_base()
      << "Using the caliper configuration: " << config << "\n";

    return config;
  }
#endif

private:
  int                      argc;
  std::vector<std::string> tokens;

  [[nodiscard]] const std::string &
  get_cmd_option(const std::string &option) const
  {
    auto itr = std::find(tokens.begin(), tokens.end(), option);
    if (itr != tokens.end())
      {
        ++itr; // Increment the iterator to get the value
        if (itr != tokens.end())
          {
            return *itr;
          }
      }
    static const std::string empty_string;
    return empty_string;
  }

  [[nodiscard]] bool
  cmd_option_exists(const std::string &option) const
  {
    return std::find(tokens.begin(), tokens.end(), option) != tokens.end();
  }
};

PRISMS_PF_END_NAMESPACE
