#ifndef parse_cmd_options_h
#define parse_cmd_options_h

#include <deal.II/base/mpi.h>

#include <algorithm>
#include <core/conditional_ostreams.h>
#include <fstream>
#include <string>
#include <vector>

class parseCMDOptions
{
public:
  parseCMDOptions(int &_argc, char **argv)
    : argc(_argc)
  {
    for (int i = 1; i < argc; ++i)
      {
        tokens.emplace_back(argv[i]);
      }
  }

  std::string
  getParametersFilename()
  {
    // Check that there aren't too many arguments
    if (argc > 3)
      {
        throw std::runtime_error(
          "Too many command line arguments. The only argument should specify "
          "the input file name. For example, `./main -i parameters.prm`");
      }

    // Default filename
    std::string parameters_filename = "parameters.prm";

    if (argc >= 2)
      {
        if (!cmdOptionExists("-i"))
          {
            throw std::runtime_error(
              "Invalid command line option. Use `-i` to specify the input file, "
              "e.g., `./main -i parameters.prm`.");
          }

        // If -i is provided, check for an argument (argc == 3)
        if (argc == 3)
          {
            parameters_filename = getCmdOption("-i");
          }
      }

    // Validate the parameters file (i.e., it ends in .prm)
    if (parameters_filename.size() < 4 ||
        parameters_filename.substr(parameters_filename.size() - 4) != ".prm")
      {
        throw std::runtime_error("The input file must have the `.prm` extension. Please "
                                 "rename the file accordingly.");
      }

    // Check if the file exists
    std::ifstream ifs_prm(parameters_filename);
    if (!ifs_prm)
      {
        throw std::runtime_error("The specified parameters file `" + parameters_filename +
                                 "` does not exist.");
      }

    // Log the filename being used
    conditionalOStreams::pout_base
      << "Using the input parameter file: " << parameters_filename << '\n';

    return parameters_filename;
  }

private:
  int                      argc;
  std::vector<std::string> tokens;

  [[nodiscard]] const std::string &
  getCmdOption(const std::string &option) const
  {
    std::vector<std::string>::const_iterator itr;
    itr = std::find(tokens.begin(), tokens.end(), option);
    if (itr != tokens.end() && ++itr != tokens.end())
      {
        return *itr;
      }
    static const std::string empty_string;
    return empty_string;
  }

  [[nodiscard]] bool
  cmdOptionExists(const std::string &option) const
  {
    return std::find(tokens.begin(), tokens.end(), option) != tokens.end();
  }
};

#endif
