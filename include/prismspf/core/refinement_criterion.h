// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef refinement_criterion_h
#define refinement_criterion_h

#include <prismspf/config.h>

#include <cfloat>
#include <cstdint>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

namespace GridRefinement
{
  /**
   * \brief Flags for refinement criterion.
   */
  enum RefinementFlags : std::uint8_t
  {
    /**
     * \brief No adaptive refinement criterion.
     */
    nothing = 0,

    /**
     * \brief Use value of the variable as a criterion for refinement.
     */
    value = 0x0001,

    /**
     * \brief Use gradient of the variable as a criterion for refinement.
     */
    gradient = 0x0002,
  };

  // Function that enables bitwise OR between flags
  inline RefinementFlags
  operator|(const RefinementFlags flag_1, const RefinementFlags flag_2)
  {
    return static_cast<RefinementFlags>(static_cast<unsigned int>(flag_1) |
                                        static_cast<unsigned int>(flag_2));
  }

  // Function that enables bitwise compound OR between flags
  inline RefinementFlags &
  operator|=(RefinementFlags &flag_1, const RefinementFlags flag_2)
  {
    flag_1 = flag_1 | flag_2;
    return flag_1;
  }

  // Function that enables bitwise AND between flags
  inline RefinementFlags
  operator&(const RefinementFlags flag_1, const RefinementFlags flag_2)
  {
    return static_cast<RefinementFlags>(static_cast<unsigned int>(flag_1) &
                                        static_cast<unsigned int>(flag_2));
  }

  // Function that enables bitwise compound AND between flags
  inline RefinementFlags &
  operator&=(RefinementFlags &flag_1, const RefinementFlags flag_2)
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
    std::string     variable_name;
    unsigned int    variable_index       = 0;
    RefinementFlags criterion_type       = RefinementFlags::nothing;
    double          value_lower_bound    = DBL_MAX;
    double          value_upper_bound    = DBL_MAX;
    double          gradient_lower_bound = DBL_MAX;

    /**
     * \brief Convert refinement criterion type to string.
     */
    [[nodiscard]] std::string
    criterion_to_string() const
    {
      if (criterion_type == RefinementFlags::nothing)
        {
          return "None";
        }
      if (((criterion_type & RefinementFlags::value) != 0U) &&
          ((criterion_type & RefinementFlags::gradient) != 0U))
        {
          return "Value and gradient";
        }
      if ((criterion_type & RefinementFlags::value) != 0U)
        {
          return "Value";
        }
      if ((criterion_type & RefinementFlags::gradient) != 0U)
        {
          return "Gradient";
        }

      return "Unknown criterion";
    }
  };
} // namespace GridRefinement

PRISMS_PF_END_NAMESPACE

#endif