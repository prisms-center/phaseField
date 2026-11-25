// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mpi.h>

#include <prismspf/utilities/assert.h>

#include <prismspf/config.h>

#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Class to parse command line options.
 */
class ParseCMDOptions
{
public:
  ParseCMDOptions(int argc, char *argv[])
  {
    int opt = 0;
    while ((opt = getopt(argc, argv, "i:o:c:h")) != -1)
      {
        switch (opt)
          {
            case 'i':
              parameters_file_name = optarg;
              break;
            case 'o':
              output_dir = optarg;
              break;
            case 'c':
              // The actual Caliper configuration options are used in the MPIInitFinalize
              // call. Here, we just reserve the option and print it in the help option.
              break;
            case 'h':
              print_usage(argv[0]);
              std::exit(EXIT_SUCCESS);
            default:
              break;
          }
      }
  };

  /**
   * @brief Print command line usage.
   */
  static void
  print_usage(const char *program_name)
  {
    std::cout << "Usage: " << program_name << " -i <input_file> -o <output_directory>\n";
    std::cout << "Options:\n";
    std::cout << "  -i <file>     Parameters file (default: parameters.prm)\n";
    std::cout << "  -o <dir>      Output directory (default: output/)\n";
#ifdef PRISMS_PF_WITH_CALIPER
    std::cout << "  -c <config>   Caliper configuration (default: "
                 "runtime-report,mem.highwatermark)\n";
#endif
    std::cout << "  -h            Show this help message\n";
  }

  [[nodiscard]] std::string
  get_parameters_filename()
  {
    // Check that the file exists and is readable
    std::ifstream file(parameters_file_name);
    ASSERT(file.is_open(),
           "The specified parameters file was unable to be opened",
           parameters_file_name);

    return parameters_file_name;
  }

private:
  std::string parameters_file_name = "parameters.prm";
  std::string output_dir           = "output/";
};

PRISMS_PF_END_NAMESPACE
