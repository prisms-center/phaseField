#ifndef refinement_criterion_h
#define refinement_criterion_h

#include <prismspf/config.h>

#include <cfloat>
#include <cstdint>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

enum RefinementCriterionFlags : std::uint8_t
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
  unsigned int             variable_index = 0;
  RefinementCriterionFlags criterion_type = RefinementCriterionFlags::criterion_default;
  double                   value_lower_bound    = DBL_MAX;
  double                   value_upper_bound    = DBL_MAX;
  double                   gradient_lower_bound = DBL_MAX;

  /**
   * \brief Convert refinement criterion type to string.
   */
  [[nodiscard]] std::string
  criterion_to_string() const
  {
    if (criterion_type == RefinementCriterionFlags::criterion_default)
      {
        return "None";
      }
    if (((criterion_type & RefinementCriterionFlags::criterion_value) != 0U) &&
        ((criterion_type & RefinementCriterionFlags::criterion_gradient) != 0U))
      {
        return "Value and gradient";
      }
    if ((criterion_type & RefinementCriterionFlags::criterion_value) != 0U)
      {
        return "Value";
      }
    if ((criterion_type & RefinementCriterionFlags::criterion_gradient) != 0U)
      {
        return "Gradient";
      }

    return "Unknown criterion";
  }
};

PRISMS_PF_END_NAMESPACE

#endif