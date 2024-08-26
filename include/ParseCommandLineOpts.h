#ifndef INCLUDE_PARSECOMMANDLINE_OPTS_H_
#define INCLUDE_PARSECOMMANDLINE_OPTS_H_

#include <deal.II/base/mpi.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

class ParseCommandLineOpts
{
public:
  ParseCommandLineOpts(int &_argc, char **argv)
  {
    argc = _argc;
    for (int i = 1; i < argc; ++i)
      tokens.push_back(std::string(argv[i]));
  }

  std::string
  getParametersFilename()
  {
    std::string parameters_filename;

    if (argc == 3)
      {
        if (cmdOptionExists("-i"))
          {
            parameters_filename = getCmdOption("-i");
            if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
              {
                std::cout << "Using the input parameter file: " << parameters_filename
                          << std::endl;
              }
          }
        else
          {
            throw("Invalid command line option given. The only argument should "
                  "be to specify the input file name.");
          }
      }
    else if (argc == 2)
      {
        if (cmdOptionExists("-i"))
          {
            parameters_filename = "parameters.prm";
            std::ifstream ifs_prm(parameters_filename);
            std::ifstream ifs_in("parameters.in");
            if (!ifs_prm && ifs_in)
              throw("The previous extension .in for the parameters file is no "
                    "longer accepted. Please rename parameters.in as "
                    "parameters.prm");
            if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
              {
                std::cout << "Using the input parameter file: " << parameters_filename
                          << std::endl;
              }
          }
        else
          {
            throw("Invalid command line option given. The only argument should "
                  "be to specify the input file name.");
          }
      }
    else if (argc == 1)
      {
        parameters_filename = "parameters.prm";
        std::ifstream ifs_prm(parameters_filename);
        std::ifstream ifs_in("parameters.in");
        if (!ifs_prm && ifs_in)
          throw("The previous extension .in for the parameters file is no longer "
                "accepted. Please rename parameters.in as parameters.prm");
        if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
          {
            std::cout << "Using the input parameter file: " << parameters_filename
                      << std::endl;
          }
      }
    else
      {
        throw("Too many command line arguments were given. The only argument "
              "should be to specify the input file name.");
      }

    return parameters_filename;
  }

private:
  int                      argc;
  std::vector<std::string> tokens;

  const std::string &
  getCmdOption(const std::string &option) const
  {
    std::vector<std::string>::const_iterator itr;
    itr = std::find(tokens.begin(), tokens.end(), option);
    if (itr != tokens.end() && ++itr != tokens.end())
      {
        return *itr;
      }
    static const std::string empty_string("");
    return empty_string;
  }

  bool
  cmdOptionExists(const std::string &option) const
  {
    return std::find(tokens.begin(), tokens.end(), option) != tokens.end();
  }
};

#endif
