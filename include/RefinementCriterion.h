#ifndef INCLUDE_REFINMENTCRITERION_H_
#define INCLUDE_REFINMENTCRITERION_H_

enum RefinementCriterionFlags
{
  criterion_default  = 0,
  criterion_value    = 0x0001,
  criterion_gradient = 0x0002
};

// Function that enables bitwise OR between flags
inline RefinementCriterionFlags
operator|(const RefinementCriterionFlags flag_1, const RefinementCriterionFlags flag_2)
{
  return static_cast<RefinementCriterionFlags>(static_cast<unsigned int>(flag_1) |
                                               static_cast<unsigned int>(flag_2));
}

// Function that enables bitwise compound OR between flags
inline RefinementCriterionFlags &
operator|=(RefinementCriterionFlags &flag_1, const RefinementCriterionFlags flag_2)
{
  flag_1 = flag_1 | flag_2;
  return flag_1;
}

// Function that enables bitwise AND between flags
inline RefinementCriterionFlags
operator&(const RefinementCriterionFlags flag_1, const RefinementCriterionFlags flag_2)
{
  return static_cast<RefinementCriterionFlags>(static_cast<unsigned int>(flag_1) &
                                               static_cast<unsigned int>(flag_2));
}

// Function that enables bitwise compound AND between flags
inline RefinementCriterionFlags &
operator&=(RefinementCriterionFlags &flag_1, const RefinementCriterionFlags flag_2)
{
  flag_1 = flag_1 & flag_2;
  return flag_1;
}

/**
 * This class holds information for a determining whether the mesh should be
 * refined.
 */
class RefinementCriterion
{
public:
  std::string              variable_name;
  unsigned int             variable_index;
  RefinementCriterionFlags criterion_type;
  double                   value_lower_bound;
  double                   value_upper_bound;
  double                   gradient_lower_bound;
};

#endif /* INCLUDE_REFINMENTCRITERION_H_ */
