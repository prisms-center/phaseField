// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>

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
    Nothing = 0,

    /**
     * \brief Use value of the variable as a criterion for refinement.
     */
    Value = 0x0001,

    /**
     * \brief Use gradient of the variable as a criterion for refinement.
     */
    Gradient = 0x0002,
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
   *
   * TODO (landinjm): This should likely be part of VariableAttributes
   */
  class RefinementCriterion
  {
  public:
    /**
     * \brief Default constructor.
     *
     * TODO (landinjm): Add some assertions here to make sure that nonsensical criteria
     * aren't specified
     */
    explicit RefinementCriterion(const RefinementFlags &_criterion,
                                 const double          &_value_lower_bound    = DBL_MAX,
                                 const double          &_value_upper_bound    = DBL_MAX,
                                 const double          &_gradient_lower_bound = DBL_MAX)
      : criterion(_criterion)
      , value_lower_bound(_value_lower_bound)
      , value_upper_bound(_value_upper_bound)
      , gradient_lower_bound(_gradient_lower_bound) {};

    /**
     * \brief Set the refinement criterion.
     */
    void
    set_criterion(const RefinementFlags &_criterion)
    {
      // TODO (landinjm): Add assertions for nontrivial getters/setters
      criterion = _criterion;
    };

    /**
     * \brief Set the value lower bound.
     */
    void
    set_value_lower_bound(const double &_value_lower_bound)
    {
      // TODO (landinjm): Add assertions for nontrivial getters/setters
      value_lower_bound = _value_lower_bound;
    };

    /**
     * \brief Set the value upper bound.
     */
    void
    set_value_upper_bound(const double &_value_upper_bound)
    {
      // TODO (landinjm): Add assertions for nontrivial getters/setters
      value_upper_bound = _value_upper_bound;
    };

    /**
     * \brief Set the gradient lower bound.
     */
    void
    set_gradient_lower_bound(const double &_gradient_lower_bound)
    {
      // TODO (landinjm): Add assertions for nontrivial getters/setters
      gradient_lower_bound = _gradient_lower_bound;
    };

    /**
     * \brief Get the refinement criterion.
     */
    [[nodiscard]] const RefinementFlags &
    get_criterion() const
    {
      // TODO (landinjm): Add assertions for nontrivial getters/setters
      return criterion;
    };

    /**
     * \brief Get the value lower bound.
     */
    [[nodiscard]] const double &
    get_value_lower_bound() const
    {
      // TODO (landinjm): Add assertions for nontrivial getters/setters
      return value_lower_bound;
    };

    /**
     * \brief Get the value upper bound.
     */
    [[nodiscard]] const double &
    get_value_upper_bound() const
    {
      // TODO (landinjm): Add assertions for nontrivial getters/setters
      return value_upper_bound;
    };

    /**
     * \brief Get the gradient lower bound.
     */
    [[nodiscard]] const double &
    get_gradient_lower_bound() const
    {
      // TODO (landinjm): Add assertions for nontrivial getters/setters
      return gradient_lower_bound;
    };

    /**
     * \brief Convert refinement criterion type to string.
     */
    [[nodiscard]] std::string
    criterion_to_string() const
    {
      if (criterion == RefinementFlags::Nothing)
        {
          return "None";
        }
      if (((criterion & RefinementFlags::Value) != 0U) &&
          ((criterion & RefinementFlags::Gradient) != 0U))
        {
          return "Value and gradient";
        }
      if ((criterion & RefinementFlags::Value) != 0U)
        {
          return "Value";
        }
      if ((criterion & RefinementFlags::Gradient) != 0U)
        {
          return "Gradient";
        }

      return "Unknown criterion";
    }

  private:
    RefinementFlags criterion            = RefinementFlags::Nothing;
    double          value_lower_bound    = DBL_MAX;
    double          value_upper_bound    = DBL_MAX;
    double          gradient_lower_bound = DBL_MAX;
  };
} // namespace GridRefinement

PRISMS_PF_END_NAMESPACE