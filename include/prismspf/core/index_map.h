
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <map>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class handles the generation and representation of each field
 * with an associated index.
 *
 * As the complexity of the types of fields we can have grows, it becomes
 * burdensome to keep track of nest maps and/or maps with pairs and tuples. This
 * class aims to unify the mapping of unique fields to some index so that we can
 * use vector to manipulate data. Additionally, it centralizes the index
 * mapping, as opposed to having local mappings spread throughout the code.
 */

class indexMap
{
public:
  /**
   *  \brief Blank constructor.
   */
  indexMap() = default;

  /**
   * \brief Default Constructor.
   */
  explicit indexMap(
    const std::map<unsigned int, variableAttributes> &_variable_attributes);

  /**
   * \brief Initialize the mapping.
   */
  void
  init();

private:
  /**
   * \brief Variable attributes.
   */
  const std::map<unsigned int, variableAttributes> *variable_attributes;

  /**
   * \brief Mapping between field index, dependency type, and the numbering in
   * indexMap.
   */
  std::vector<std::pair<unsigned int, dependencyType>> mapping;
};

indexMap::indexMap(const std::map<unsigned int, variableAttributes> &_variable_attributes)
  : variable_attributes(&_variable_attributes)
{
  init();
}

void
indexMap::init()
{
  // Loop through the variable_attributes and fill in the mappings. This is done
  // in such a way that the vector is ordered in terms of index and then
  // dependencyType. For example, 1 NORMAL, 1 OLD_1, 2 NORMAL, 3 NORMAL, 3
  // OLD_1, 3 OLD_2...

  // For now, I'll do a double loop through the variableAttributes to get the
  // right ordering.
  for (const auto &[index_1, variable_1] : *variable_attributes)
    {
      mapping.push_back(std::make_pair(index_1, dependencyType::NORMAL));

      for (const auto &[index_2, variable_2] : *variable_attributes)
        {
        }
    }
}

PRISMS_PF_END_NAMESPACE
