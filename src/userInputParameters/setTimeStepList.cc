#include <deal.II/base/exceptions.h>

#include "../../include/userInputParameters.h"
#include <vector>

template <int dim>
std::vector<unsigned int>
userInputParameters<dim>::setTimeStepList(
  const std::string               &outputSpacingType,
  unsigned int                     numberOfOutputs,
  const std::vector<unsigned int> &userGivenTimeStepList)
{
  // Initialize timestep list
  std::vector<unsigned int> timeStepList;

  // The number of outputs cannot be greater than the number increments
  if (numberOfOutputs > totalIncrements)
    {
      numberOfOutputs = totalIncrements;
    }

  // Prevent divide by zero in subsequent output types by returning the a vector where the
  // only entry is one greater than the number of increments. This way, we effectively
  // have no outputs. While this condition can be ignored for the LIST type, the user
  // should just ignore the parameter `set Number of outputs` and use the default value
  // of 10.
  if (numberOfOutputs == 0)
    {
      timeStepList.push_back(totalIncrements + 1);
      return timeStepList;
    }

  // Set output list for all the output list types
  if (outputSpacingType == "LIST")
    {
      timeStepList = userGivenTimeStepList;
    }
  else if (outputSpacingType == "EQUAL_SPACING")
    {
      for (unsigned int iter = 0; iter <= totalIncrements;
           iter += totalIncrements / numberOfOutputs)
        {
          timeStepList.push_back(iter);
        }
    }
  else if (outputSpacingType == "LOG_SPACING")
    {
      timeStepList.push_back(0);
      for (unsigned int output = 1; output <= numberOfOutputs; output++)
        {
          timeStepList.push_back(round(std::pow(static_cast<double>(totalIncrements),
                                                static_cast<double>(output) /
                                                  static_cast<double>(numberOfOutputs))));
        }
    }
  else if (outputSpacingType == "N_PER_DECADE")
    {
      AssertThrow(totalIncrements > 1,
                  dealii::ExcMessage(
                    std::string("PRISMS-PF Error: For n per decaded spaced outputs, "
                                "the number of increments must be greater than 1.")));

      timeStepList.push_back(0);
      timeStepList.push_back(1);
      for (unsigned int iter = 2; iter <= totalIncrements; iter++)
        {
          unsigned int decade    = std::ceil(std::log10(iter));
          unsigned int step_size = std::pow(10, decade) / numberOfOutputs;
          if (iter % step_size == 0)
            {
              timeStepList.push_back(iter);
            }
        }
    }
  else
    {
      AssertThrow(false,
                  dealii::ExcMessage(
                    std::string("PRISMS-PF Error: Invalid output spacing type.")));
    }

  return timeStepList;
}

// Template instantiations
#include "../../include/userInputParameters_template_instantiations.h"
