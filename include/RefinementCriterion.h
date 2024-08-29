#ifndef INCLUDE_REFINMENTCRITERION_H_
#define INCLUDE_REFINMENTCRITERION_H_

enum criterionType
{
  VALUE,
  GRADIENT,
  VALUE_AND_GRADIENT
};

/**
 * This class holds information for a determining whether the mesh should be
 * refined.
 */
class RefinementCriterion
{
public:
  std::string   variable_name;
  unsigned int  variable_index;
  criterionType criterion_type;
  double        value_lower_bound;
  double        value_upper_bound;
  double        gradient_lower_bound;
};

#endif /* INCLUDE_REFINMENTCRITERION_H_ */
