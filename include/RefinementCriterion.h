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
operator|(const RefinementCriterionFlags f1, const RefinementCriterionFlags f2)
{
  return static_cast<RefinementCriterionFlags>(static_cast<unsigned int>(f1) |
                                               static_cast<unsigned int>(f2));
}

// Function that enables bitwise compound OR between flags
inline RefinementCriterionFlags &
operator|=(RefinementCriterionFlags &f1, const RefinementCriterionFlags f2)
{
  f1 = f1 | f2;
  return f1;
}

// Function that enables bitwise AND between flags
inline RefinementCriterionFlags
operator&(const RefinementCriterionFlags f1, const RefinementCriterionFlags f2)
{
  return static_cast<RefinementCriterionFlags>(static_cast<unsigned int>(f1) &
                                               static_cast<unsigned int>(f2));
}

// Function that enables bitwise compound AND between flags
inline RefinementCriterionFlags &
operator&=(RefinementCriterionFlags &f1, const RefinementCriterionFlags f2)
{
  f1 = f1 & f2;
  return f1;
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
