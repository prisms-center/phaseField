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

  static unsigned int
  calculate_number_of_events(const double &rate,
                             const double &delta_t,
                             const double &volume,
                             RNGEngine    &rng)
  {
    std::poisson_distribution<unsigned int> distribution(rate * delta_t * volume);
    return distribution(rng);
  }

  static void
  mpi_gather_nuclei(std::vector<Nucleus<dim>> &local_nuclei);

  static void
  mpi_broadcast_nuclei(std::vector<Nucleus<dim>> &local_nuclei);
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
          std::list<Nucleus<dim>> new_nuclei_list;
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

                      new_nuclei_list.emplace_back(nucleating_index,
                                                   nucleus_location,
                                                   seed_time,
                                                   seed_increment);
                    }
                }
            }
          if (dealii::Utilities::MPI::sum(new_nuclei_list.size(),
                                          MPI_COMM_WORLD)) // dont waste time if
                                                           // no nuclei appeared
            {
              // Gather to root process
              std::vector<Nucleus<dim>> new_nuclei(new_nuclei_list.begin(),
                                                   new_nuclei_list.end());
              mpi_gather_nuclei(new_nuclei);
              if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
                {
                  // Remove nuclei within their exclusion distance and add to nuclei list
                  ConditionalOStreams::pout_base()
                    << new_nuclei.size() << " nuclei generated before exclusion.\n"
                    << "Excluding nuclei...\n";
                  unsigned int count = 0;

                  // remove bias from cell order
                  std::shuffle(new_nuclei.begin(), new_nuclei.end(), rng);

                  while (!new_nuclei.empty())
                    {
                      Nucleus<dim> *nuc = &new_nuclei.back();
                      for (const auto &existing_nucleus : nuclei)
                        {
                          if (nuc->location.distance(existing_nucleus.location) <
                                nucleation_parameters.get_exclusion_distance() ||
                              (nuc->field_index == existing_nucleus.field_index &&
                               nuc->location.distance(existing_nucleus.location) <
                                 nucleation_parameters
                                   .get_same_field_exclusion_distance()))
                            {
                              new_nuclei.pop_back();
                              nuc = nullptr;
                              break;
                            }
                        }
                      if (nuc != nullptr)
                        {
                          // Note: Using push_back() in a loop is not good use for
                          // vectors. We also don't want to use reserve() on the highest
                          // possible amount because that could allocate much more space
                          // than needed. I originally was using a std::list to avoid this
                          // issue, but that is unfriendly to the MPI functions.
                          // One solution could be to convert between data structures
                          // as needed, but that also adds overhead. For now, I will
                          // assume that the total number of added nuclei is not enough to
                          // cause significant performance issues.
                          nuclei.push_back(*nuc);
                          new_nuclei.pop_back();
                          ++count;
                        }
                    }
                  ConditionalOStreams::pout_base()
                    << count << " nuclei generated after exclusion.\n"
                    << nuclei.size() << " total nuclei.\n";
                }
              mpi_broadcast_nuclei(nuclei);
            }
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
NucleationHandler<dim, degree, number>::mpi_gather_nuclei(
  std::vector<Nucleus<dim>> &local_nuclei)
{
  // Step 1: Share how many nuclei each rank has
  int              rank        = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  int              num_procs   = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  int              local_count = local_nuclei.size();
  std::vector<int> nuclei_counts_per_rank;
  if (rank == 0)
    nuclei_counts_per_rank.resize(num_procs);

  MPI_Gather(&local_count,
             1,
             MPI_INT,
             nuclei_counts_per_rank.data(),
             1,
             MPI_INT,
             0,
             MPI_COMM_WORLD);

  // Step 2: Compute displacements and allocate receive buffer on root
  std::vector<int>          recv_displacements;
  std::vector<Nucleus<dim>> gathered_nuclei;
  if (rank == 0)
    {
      recv_displacements.resize(num_procs);
      recv_displacements[0] = 0;
      for (int r = 1; r < num_procs; ++r)
        recv_displacements[r] = recv_displacements[r - 1] + nuclei_counts_per_rank[r - 1];

      int total_count = recv_displacements.back() + nuclei_counts_per_rank.back();
      gathered_nuclei.resize(total_count);
    }

  // Step 3: Gather all nuclei into root's buffer
  MPI_Gatherv(local_nuclei.data(),
              local_count,
              Nucleus<dim>::mpi_datatype(),
              gathered_nuclei.data(),
              nuclei_counts_per_rank.data(),
              recv_displacements.data(),
              Nucleus<dim>::mpi_datatype(),
              0,
              MPI_COMM_WORLD);
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
NucleationHandler<dim, degree, number>::mpi_broadcast_nuclei(
  std::vector<Nucleus<dim>> &local_nuclei)
{
  int rank  = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  int count = local_nuclei.size();
  MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0)
    local_nuclei.resize(count);

  MPI_Bcast(local_nuclei.data(), count, Nucleus<dim>::mpi_datatype(), 0, MPI_COMM_WORLD);
}

PRISMS_PF_END_NAMESPACE