// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/types.h>

#include <prismspf/user_inputs/nucleation_parameters.h>

#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>
#include <prismspf/nucleation/nucleus.h>

#include <list>
#include <map>
#include <set>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief The class handles the stochastic nucleation
 * in PRISMS-PF.
 */
template <unsigned int dim, unsigned int degree, typename number>
class NucleationHandler
{
public:
  /**
   * @brief Construct a new Nucleation Handler object
   */
  NucleationHandler();

  /**
   * @brief Destroy the Nucleation Handler object
   */
  ~NucleationHandler();

  void
  attempt_nucleation(const SolverContext<dim, degree, number> &solver_context,
                     std::vector<Nucleus<dim>>                &nuclei);

private:
  // Add private member variables and functions as needed
};

inline void
NucleationHandler::attempt_nucleation(
  const SolverContext<dim, degree, number> &solver_context,
  std::vector<Nucleus<dim>>                &nuclei)
{
  NucleationParameters &nucleation_parameters =
    solver_context.user_inputs->get_nucleation_parameters();
  for (const auto &[index, variable] :
       solver_context.user_inputs->get_variable_attributes())
    {
      if (variable.is_nucleation_rate)
        {
          std::list<Nucleus<dim>> new_nuclei;
          // Perform nucleation logic here
          // This is where you would check conditions and create nuclei
          // for (const auto& cell)

          dealii::Point<dim> nucleus_location; // Determine the location based on your
                                               // logic
          double seed_time =
            solver_context.user_inputs->get_temporal_discretization().get_time();
          unsigned int seed_increment =
            solver_context.user_inputs->get_temporal_discretization().get_increment();

          nuclei.emplace_back(index, nucleus_location, seed_time, seed_increment);
        }
    }
}

PRISMS_PF_END_NAMESPACE