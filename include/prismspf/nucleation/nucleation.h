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

  static double
  calculate_nucleation_probability(const double &nucleation_rate, const double &delta_t)
  {
    return 1.0 - std::exp(-nucleation_rate * delta_t);
  }

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
  double delta_t =
    nucleation_parameters.get_nucleation_period() *
    solver_context.user_inputs->get_temporal_discretization().get_delta_t();
  auto &rng = solver_context.user_inputs->get_miscellaneous_parameters().rng;
  for (const auto &[index, variable] :
       solver_context.user_inputs->get_variable_attributes())
    {
      if (variable.is_nucleation_rate)
        {
          std::list<Nucleus<dim>> new_nuclei;
          // Perform nucleation logic here
          // This is where you would check conditions and create nuclei
          for (const auto &cell : solver_context.get_triangulation_handler()
                                    .get_triangulation()
                                    .active_cell_iterators())
            {
              if (cell->is_locally_owned())
                {
                  // Grab the DoFHandler iterator
                  const auto dof_iterator = cell->as_dof_handler_iterator(
                    solver_context.get_dof_handler().get_dof_handler(index));

                  // Reinit the cell
                  fe_values.reinit(dof_iterator);
                  // Get the values for a scalar field
                  fe_values.get_function_values(
                    solver_context.get_solution_handler()
                      .get_solution_vector(index, DependencyType::Normal),
                    values);
                  for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                    {
                      nuc_rate +=
                        values[q_point] * fe_values.get_quadrature().weight(q_point);
                    }
                  double nuc_probability =
                    nucleation_parameters.calculate_nucleation_probability(nuc_rate,
                                                                           delta_t);
                  if (random() < nuc_probability)
                    {
                      dealii::Point<dim> nucleus_location; // mapping(Random)
                      double             seed_time =
                        solver_context.user_inputs->get_temporal_discretization()
                          .get_time();
                      unsigned int seed_increment =
                        solver_context.user_inputs->get_temporal_discretization()
                          .get_increment();
                      unsigned int nucleating_index =
                        0; // pick randomly from allowed indices
                      if (randombool)
                        {
                          new_nuclei.emplace_back(nucleating_index,
                                                  nucleus_location,
                                                  seed_time,
                                                  seed_increment);
                        }
                      else
                        {
                          new_nuclei.emplace_front(nucleating_index,
                                                   nucleus_location,
                                                   seed_time,
                                                   seed_increment);
                        }
                    }
                }
            }
        }
    }
}

PRISMS_PF_END_NAMESPACE