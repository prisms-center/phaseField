// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/fe/fe_values.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/nucleation_parameters.h>

#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>
#include <prismspf/nucleation/nucleus.h>

#include <algorithm>
#include <list>
#include <map>
#include <random>
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
  static void
  attempt_nucleation(const SolverContext<dim, degree, number> &solver_context,
                     std::vector<Nucleus<dim>>                &nuclei);

  static double
  calculate_number_of_events(const double &rate,
                             const double &delta_t,
                             const double &volume,
                             RNGEngine    &rng)
  {
    std::poisson_distribution<unsigned int> distribution(rate * delta_t * volume);
    return distribution(rng);
  }
};

template <unsigned int dim, unsigned int degree, typename number>
inline void
NucleationHandler<dim, degree, number>::attempt_nucleation(
  const SolverContext<dim, degree, number> &solver_context,
  std::vector<Nucleus<dim>>                &nuclei)
{
  const NucleationParameters &nucleation_parameters =
    solver_context.get_user_inputs().get_nucleation_parameters();
  double delta_t =
    nucleation_parameters.get_nucleation_period() *
    solver_context.get_user_inputs().get_temporal_discretization().get_timestep();
  auto &rng = solver_context.get_user_inputs().get_miscellaneous_parameters().rng;
  // Set up FEValues
  const dealii::QGaussLobatto<dim> quadrature(degree + 1);
  const unsigned int               num_quad_points = quadrature.size();
  dealii::FESystem<dim> fe_system(dealii::FE_Q<dim>(dealii::QGaussLobatto<1>(degree + 1)),
                                  1);
  dealii::FEValues<dim> fe_values(fe_system,
                                  quadrature,
                                  dealii::UpdateFlags::update_values |
                                    dealii::UpdateFlags::update_JxW_values);
  for (const auto &[index, variable] :
       solver_context.get_user_inputs().get_variable_attributes())
    {
      if (variable.is_nucleation_rate())
        {
          std::uniform_int_distribution<unsigned int> nucleating_index_dist(
            0,
            variable.get_nucleating_field_indices().size() - 1);
          std::list<Nucleus<dim>> new_nuclei;
          // Perform nucleation logic here
          // This is where you would check conditions and create nuclei
          for (const auto &cell : solver_context.get_triangulation_handler()
                                    .get_triangulation()
                                    .active_cell_iterators())
            {
              if (cell->is_locally_owned())
                {
                  std::vector<number> values(num_quad_points, 0.0);
                  // Grab the DoFHandler iterator
                  const auto dof_iterator = cell->as_dof_handler_iterator(
                    solver_context.get_dof_handler().get_dof_handler(index));

                  // Reinit the cell
                  fe_values.reinit(dof_iterator);
                  // Get the values for a scalar field
                  fe_values.get_function_values(
                    *(solver_context.get_solution_handler()
                        .get_solution_vector(index, DependencyType::Normal)),
                    values);
                  double nuc_rate    = 0.0;
                  double cell_volume = 0.0;
                  for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                    {
                      nuc_rate +=
                        values[q_point] * fe_values.get_quadrature().weight(q_point);
                      cell_volume += fe_values.JxW(q_point);
                    }
                  unsigned int num_nuclei_in_cell = calculate_number_of_events(
                    nuc_rate,
                    delta_t,
                    cell_volume,
                    solver_context.get_user_inputs().get_miscellaneous_parameters().rng);
                  for (unsigned int i = 0; i < num_nuclei_in_cell; ++i)
                    {
                      dealii::Point<dim> nucleus_location_unit_cell;
                      for (unsigned int d = 0; d < dim; ++d)
                        {
                          static std::uniform_real_distribution<double>
                            uniform_unit_interval(-1.0, 1.0);
                          nucleus_location_unit_cell[d] = uniform_unit_interval(rng);
                        }
                      dealii::Point<dim> nucleus_location =
                        solver_context.get_mapping()
                          .transform_unit_to_real_cell(cell, nucleus_location_unit_cell);
                      double seed_time = solver_context.get_user_inputs()
                                           .get_temporal_discretization()
                                           .get_time();
                      unsigned int seed_increment = solver_context.get_user_inputs()
                                                      .get_temporal_discretization()
                                                      .get_increment();
                      unsigned int nucleating_index =
                        variable
                          .get_nucleating_field_indices()[nucleating_index_dist(rng)];

                      new_nuclei.emplace_back(nucleating_index,
                                              nucleus_location,
                                              seed_time,
                                              seed_increment);
                    }
                }
            }
          // Remove nuclei within their exclusion distance and add to nuclei list
          ConditionalOStreams::pout_base()
            << new_nuclei.size() << " nuclei generated before exclusion.\n"
            << "Excluding nuclei...\n";

          // remove bias
          std::vector<Nucleus<dim>> vec(new_nuclei.begin(), new_nuclei.end());
          std::shuffle(vec.begin(), vec.end(), rng);
          new_nuclei = std::list<Nucleus<dim>>(vec.begin(), vec.end());

          while (!new_nuclei.empty())
            {
              auto iter = new_nuclei.begin();
              for (const auto &existing_nucleus : nuclei)
                {
                  if (iter->location.distance(existing_nucleus.location) <
                        nucleation_parameters.get_exclusion_distance() ||
                      (iter->field_index == existing_nucleus.field_index &&
                       iter->location.distance(existing_nucleus.location) <
                         nucleation_parameters.get_same_field_exclusion_distance()))
                    {
                      new_nuclei.erase(iter);
                      iter = new_nuclei.end();
                      break;
                    }
                }
              if (iter != new_nuclei.end())
                {
                  nuclei.push_back(*iter);
                  new_nuclei.erase(iter);
                }
            }
          ConditionalOStreams::pout_base()
            << new_nuclei.size() << " nuclei generated after exclusion.\n";
        }
    }
}

PRISMS_PF_END_NAMESPACE