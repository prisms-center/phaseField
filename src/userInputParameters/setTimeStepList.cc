#include "../../include/userInputParameters.h"

template <int dim>
std::vector<unsigned int>
userInputParameters<dim>::setTimeStepList(
  const std::string                outputSpacingType,
  unsigned int                     numberOfOutputs,
  const std::vector<unsigned int> &userGivenTimeStepList)
{
  std::vector<unsigned int> timeStepList;

  if (numberOfOutputs > 0)
    {
      if (outputSpacingType == "EQUAL_SPACING")
        {
          if (numberOfOutputs > totalIncrements)
            numberOfOutputs = totalIncrements;

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
              timeStepList.push_back(
                round(std::pow(10,
                               double(output) / double(numberOfOutputs) *
                                 std::log10(totalIncrements))));
            }
        }
      else if (outputSpacingType == "N_PER_DECADE")
        {
          timeStepList.push_back(0);
          timeStepList.push_back(1);
          for (unsigned int iter = 2; iter <= totalIncrements; iter++)
            {
              int decade    = std::ceil(std::log10(iter));
              int step_size = (std::pow(10, decade)) / numberOfOutputs;
              if (iter % step_size == 0)
                {
                  timeStepList.push_back(iter);
                }
            }
        }
    }
  else if (outputSpacingType == "LIST")
    {
      timeStepList = userGivenTimeStepList;
    }
  else
    {
      // I'm not sure why this is set up this way. It seems like the intuitive
      // thing would be to have an empty list, not a list with one entry that is
      // higher than will be reached during time stepping.
      timeStepList.push_back(totalIncrements + 1);
    }

  return timeStepList;
}

// Template instantiations
#include "../../include/userInputParameters_template_instantiations.h"
