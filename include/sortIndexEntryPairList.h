// At some point, rewrite the following four methods as a single templated
// method
#ifndef INCLUDE_SORTINDEXENTRYPAIRLIST_H_
#define INCLUDE_SORTINDEXENTRYPAIRLIST_H_

#include "varTypeEnums.h"

#include <string>
#include <vector>

std::vector<bool>
sortIndexEntryPairList(
  const std::vector<std::pair<unsigned int, bool>> &unsorted_pair_list,
  unsigned int                                      number_of_variables,
  bool                                              default_value);
std::vector<std::string>
sortIndexEntryPairList(
  const std::vector<std::pair<unsigned int, std::string>> &unsorted_pair_list,
  unsigned int                                             number_of_variables,
  const std::string                                       &default_value);
std::vector<fieldType>
sortIndexEntryPairList(
  const std::vector<std::pair<unsigned int, fieldType>> &unsorted_pair_list,
  unsigned int                                           number_of_variables,
  fieldType                                              default_value);
std::vector<PDEType>
sortIndexEntryPairList(
  const std::vector<std::pair<unsigned int, PDEType>> &unsorted_pair_list,
  unsigned int                                         number_of_variables,
  PDEType                                              default_value);

#endif
