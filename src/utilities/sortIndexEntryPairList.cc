#include "../../include/sortIndexEntryPairList.h"

std::vector<bool>
sortIndexEntryPairList(
  const std::vector<std::pair<unsigned int, bool>> &unsorted_pair_list,
  unsigned int                                      number_of_variables,
  bool                                              default_value)
{
  std::vector<bool> sorted_vec;
  unsigned int      entry_index;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      bool entry_found = false;
      for (unsigned int j = 0; j < unsorted_pair_list.size(); j++)
        {
          if (i == unsorted_pair_list.at(j).first)
            {
              entry_found = true;
              entry_index = j;
              break;
            }
        }

      if (entry_found)
        {
          sorted_vec.push_back(unsorted_pair_list.at(entry_index).second);
        }
      else
        {
          sorted_vec.push_back(default_value);
        }
    }
  return sorted_vec;
}

std::vector<std::string>
sortIndexEntryPairList(
  const std::vector<std::pair<unsigned int, std::string>> &unsorted_pair_list,
  unsigned int                                             number_of_variables,
  const std::string                                       &default_value)
{
  std::vector<std::string> sorted_vec;
  unsigned int             entry_index;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      bool entry_found = false;
      for (unsigned int j = 0; j < unsorted_pair_list.size(); j++)
        {
          if (i == unsorted_pair_list.at(j).first)
            {
              entry_found = true;
              entry_index = j;
              break;
            }
        }

      if (entry_found)
        {
          sorted_vec.push_back(unsorted_pair_list.at(entry_index).second);
        }
      else
        {
          sorted_vec.push_back(default_value);
        }
    }
  return sorted_vec;
}

std::vector<fieldType>
sortIndexEntryPairList(
  const std::vector<std::pair<unsigned int, fieldType>> &unsorted_pair_list,
  unsigned int                                           number_of_variables,
  fieldType                                              default_value)
{
  std::vector<fieldType> sorted_vec;
  unsigned int           entry_index;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      bool entry_found = false;
      for (unsigned int j = 0; j < unsorted_pair_list.size(); j++)
        {
          if (i == unsorted_pair_list.at(j).first)
            {
              entry_found = true;
              entry_index = j;
              break;
            }
        }

      if (entry_found)
        {
          sorted_vec.push_back(unsorted_pair_list.at(entry_index).second);
        }
      else
        {
          sorted_vec.push_back(default_value);
        }
    }
  return sorted_vec;
}

std::vector<PDEType>
sortIndexEntryPairList(
  const std::vector<std::pair<unsigned int, PDEType>> &unsorted_pair_list,
  unsigned int                                         number_of_variables,
  PDEType                                              default_value)
{
  std::vector<PDEType> sorted_vec;
  unsigned int         entry_index;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      bool entry_found = false;
      for (unsigned int j = 0; j < unsorted_pair_list.size(); j++)
        {
          if (i == unsorted_pair_list.at(j).first)
            {
              entry_found = true;
              entry_index = j;
              break;
            }
        }

      if (entry_found)
        {
          sorted_vec.push_back(unsorted_pair_list.at(entry_index).second);
        }
      else
        {
          sorted_vec.push_back(default_value);
        }
    }
  return sorted_vec;
}
