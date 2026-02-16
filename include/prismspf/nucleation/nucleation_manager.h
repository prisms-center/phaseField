// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mpi.h>
#include <deal.II/fe/fe_values.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/phase_field_tools.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/miscellaneous_parameters.h>
#include <prismspf/user_inputs/nucleation_parameters.h>
#include <prismspf/user_inputs/temporal_discretization.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/solve_context.h>

#include <prismspf/utilities/periodic_distance.h>

#include <prismspf/config.h>
#include <prismspf/nucleation/nucleus.h>

#include <algorithm>
#include <list>
#include <mpi.h>
#include <random>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

// Note: could make NucleationManager a namespace as everything is static

/**
 * @brief The class handles the stochastic nucleation
 * in PRISMS-PF.
 */
template <unsigned int dim, unsigned int degree, typename number>
class NucleationManager
{
public:
  /**
   * @brief Main nucleation function. Iterates over the domain and stochastically adds
   * nuclei to the list @param nuclei based on nucleation rates.
   */
  static bool
  attempt_nucleation(const SolveContext<dim, degree, number> &solve_context,
                     std::vector<Nucleus<dim>>               &nuclei);

  /**
   * @brief Samples the poisson distribution to calculate a number of events in a
   * time-volume given a rate.
   */
  static unsigned int
  calculate_number_of_events(const double &rate,
                             const double &delta_t,
                             const double &volume,
                             RNGEngine    &rng)
  {
    std::poisson_distribution<unsigned int> distribution(rate * delta_t * volume);
    return distribution(rng);
  }

  /**
   * @brief Gathers the potential new nuclei from each processor onto the root process,
   * eliminates any nuclei that need be excluded, then updates and redistributes the
   * global list to each process.
   */
  static bool
  gather_exclude_broadcast_nuclei(std::list<Nucleus<dim>>        &new_nuclei_list,
                                  std::vector<Nucleus<dim>>      &global_nuclei,
                                  const UserInputParameters<dim> &user_inputs);

  /**
   * @brief Gathers nuclei lists to root.
   * Modifies @param local_nuclei
   */
  static void
  mpi_gather_nuclei(std::vector<Nucleus<dim>> &local_nuclei);

  /**
   * @brief Broadcasts nuclei lists from root.
   * Modifies @param local_nuclei
   */
  static void
  mpi_broadcast_nuclei(std::vector<Nucleus<dim>> &local_nuclei);
};

template <unsigned int dim, unsigned int degree, typename number>
inline bool
NucleationManager<dim, degree, number>::attempt_nucleation(
  const SolveContext<dim, degree, number> &solve_context,
  std::vector<Nucleus<dim>>               &nuclei)
{
  // Set up references.
  const UserInputParameters<dim> &user_inputs = solve_context.get_user_inputs();
  const NucleationParameters     &nuc_params  = user_inputs.get_nucleation_parameters();
  const TemporalDiscretization   &time_info   = user_inputs.get_temporal_discretization();
  const double delta_t = nuc_params.get_nucleation_period() * time_info.get_timestep();
  auto        &rng     = user_inputs.get_miscellaneous_parameters().rng;

  // Set up FEValues
  static const dealii::QGaussLobatto<dim> quadrature(degree + 1);
  static const unsigned int               num_quad_points = quadrature.size();
  dealii::FESystem<dim> fe_system(dealii::FE_Q<dim>(dealii::QGaussLobatto<1>(degree + 1)),
                                  1);
  dealii::FEValues<dim> fe_values(fe_system,
                                  quadrature,
                                  dealii::UpdateFlags::update_values |
                                    dealii::UpdateFlags::update_JxW_values);
  std::list<Nucleus<dim>> new_nuclei_list;
  // Loop over nucleation rate variables and attempt seeding at each cell
  for (unsigned int index = 0; index < solve_context.get_field_attributes().size();
       ++index)
    {
      const auto &variable = solve_context.get_field_attributes()[index];
      if (!variable.is_nucleation_rate_variable)
        {
          continue;
        }
      std::uniform_int_distribution<unsigned int> nucleating_index_dist(
        0,
        variable.nucleating_field_indices.size() - 1);
      // Perform nucleation logic here
      // This is where you would check conditions and create nuclei
      for (const auto &cell : solve_context.get_triangulation_manager()
                                .get_triangulation()
                                .active_cell_iterators())
        {
          if (!cell->is_locally_owned())
            {
              continue;
            }
          std::vector<number> values(num_quad_points, 0.0);
          // Grab the DoFHandler iterator
          const auto dof_iterator = cell->as_dof_handler_iterator(
            solve_context.get_dof_manager().get_field_dof_handler(index));

          // Reinit the cell
          fe_values.reinit(dof_iterator);
          // Get the values for a scalar field
          fe_values.get_function_values(
            (solve_context.get_solution_indexer().get_solution_vector(index)),
            values);
          double nuc_rate    = 0.0;
          double cell_volume = 0.0;
          for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
            {
              nuc_rate += values[q_point] * fe_values.get_quadrature().weight(q_point);
              cell_volume += fe_values.JxW(q_point);
            }
          unsigned int num_nuclei_in_cell =
            calculate_number_of_events(nuc_rate, delta_t, cell_volume, rng);
          for (unsigned int i = 0; i < num_nuclei_in_cell; ++i)
            {
              dealii::Point<dim> nucleus_location_unit_cell;
              for (unsigned int d = 0; d < dim; ++d)
                {
                  // Note: if we ever do non-rectangular cells, try a randomly weighted
                  // sum over dealii::GeometryInfo< dim >::unit_cell_vertex
                  static std::uniform_real_distribution<double> uniform_unit_interval(
                    0.0,
                    1.0);
                  nucleus_location_unit_cell[d] = uniform_unit_interval(rng);
                }
              dealii::Point<dim> nucleus_location =
                SystemWide<dim, 0>::mapping
                  .transform_unit_to_real_cell(cell, nucleus_location_unit_cell);
              double       seed_time      = time_info.get_time();
              unsigned int seed_increment = time_info.get_increment();
              unsigned int nucleating_index =
                variable.nucleating_field_indices[nucleating_index_dist(rng)];

              new_nuclei_list.emplace_back(nucleating_index,
                                           nucleus_location,
                                           seed_time,
                                           seed_increment);
            }
        }
    }
  return gather_exclude_broadcast_nuclei(new_nuclei_list, nuclei, user_inputs);
  //// If any nuclei were created, gather them to rank 0, exclude based on distance, and
  //// broadcast back to all ranks
  // if (dealii::Utilities::MPI::sum(new_nuclei_list.size(),
  //                                 MPI_COMM_WORLD)) // dont waste time if
  //                                                  // no nuclei appeared
  //   {
  //     // Gather to root process
  //     std::vector<Nucleus<dim>> new_nuclei(new_nuclei_list.begin(),
  //                                          new_nuclei_list.end());
  //     mpi_gather_nuclei(new_nuclei);
  //     if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  //       {
  //         // Remove nuclei within their exclusion distance and add to nuclei list
  //         ConditionalOStreams::pout_base()
  //           << new_nuclei.size() << " nuclei generated before exclusion.\n"
  //           << "Excluding nuclei...\n";
  //         unsigned int count = 0;
  //
  //        // remove bias from cell order
  //        std::shuffle(new_nuclei.begin(), new_nuclei.end(), rng);
  //
  //        while (!new_nuclei.empty())
  //          {
  //            Nucleus<dim> &nuc   = new_nuclei.back();
  //            bool          valid = std::none_of(
  //              nuclei.begin(),
  //              nuclei.end(),
  //              [&](const Nucleus<dim> &existing_nucleus)
  //              {
  //                const double distance =
  //                  prisms::distance<dim, double>(nuc.location,
  //                                                existing_nucleus.location,
  //                                                user_inputs);
  //                return nuc_params.check_active(existing_nucleus, time_info) &&
  //                       (distance < nuc_params.get_exclusion_distance() ||
  //                        (nuc.field_index == existing_nucleus.field_index &&
  //                         distance < nuc_params.get_same_field_exclusion_distance()));
  //              });
  //            if (valid)
  //              {
  //                // Note: Using push_back() in a loop is not good use for
  //                // vectors. We also don't want to use reserve() on the upper bound
  //                // because that could allocate much more space than needed. I
  //                originally
  //                // was using a std::list to avoid this issue, but that is unfriendly
  //                to
  //                // the MPI functions. One solution could be to convert between data
  //                // structures as needed, but that also adds overhead. For now, I will
  //                // assume that the total number of nuclei is not enough to
  //                // cause significant performance issues.
  //                nuclei.push_back(nuc);
  //                ++count;
  //                any_nucleation_occurred = true;
  //              }
  //            new_nuclei.pop_back();
  //          }
  //        ConditionalOStreams::pout_base()
  //          << count << " nuclei generated after exclusion.\n"
  //          << nuclei.size() << " total nuclei.\n";
  //      }
  //    MPI_Bcast(&any_nucleation_occurred, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
  //    mpi_broadcast_nuclei(nuclei);
  //  }
  // return any_nucleation_occurred;
}

template <unsigned int dim, unsigned int degree, typename number>
inline bool
NucleationManager<dim, degree, number>::gather_exclude_broadcast_nuclei(
  std::list<Nucleus<dim>>        &new_nuclei_list,
  std::vector<Nucleus<dim>>      &global_nuclei,
  const UserInputParameters<dim> &user_inputs)
{
  // dont waste time if no nuclei appeared
  if (!bool(dealii::Utilities::MPI::sum(new_nuclei_list.size(), MPI_COMM_WORLD)))
    {
      return false;
    }

  // Set up refs
  const NucleationParameters   &nuc_params = user_inputs.get_nucleation_parameters();
  const TemporalDiscretization &time_info  = user_inputs.get_temporal_discretization();
  RNGEngine                    &rng = user_inputs.get_miscellaneous_parameters().rng;

  // Gather new nuclei to root process
  std::vector<Nucleus<dim>> new_nuclei(new_nuclei_list.begin(), new_nuclei_list.end());
  mpi_gather_nuclei(new_nuclei);
  bool any_nucleation_occurred = false;
  if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      // Remove nuclei within their exclusion distance and add to nuclei list
      ConditionalOStreams::pout_base()
        << "[Increment " << time_info.get_increment() << "] : Nucleation\n"
        << "  " << new_nuclei.size() << " nuclei generated before exclusion.\n"
        << "  Excluding nuclei...\n";
      unsigned int count = 0;

      // remove bias from cell order
      std::shuffle(new_nuclei.begin(), new_nuclei.end(), rng);

      while (!new_nuclei.empty())
        {
          const Nucleus<dim> &nuc   = new_nuclei.back();
          bool                valid = std::none_of(
            global_nuclei.begin(),
            global_nuclei.end(),
            [&](const Nucleus<dim> &existing_nucleus)
              {
                const double distance =
                  prisms::distance<dim, double>(nuc.location,
                                                existing_nucleus.location,
                                                user_inputs);
                return nuc_params.check_active(existing_nucleus, time_info) &&
                       (distance < nuc_params.get_exclusion_distance() ||
                        (nuc.field_index == existing_nucleus.field_index &&
                         distance < nuc_params.get_same_field_exclusion_distance()));
              });
          if (valid)
            {
              // Note: Using push_back() in a loop is not good use for
              // vectors. We also don't want to use reserve() on the upper bound
              // because that could allocate much more space than needed. I originally
              // was using a std::list to avoid this issue, but that is unfriendly to
              // the MPI functions. One solution could be to convert between data
              // structures as needed, but that also adds overhead. For now, I will
              // assume that the total number of nuclei is not enough to
              // cause significant performance issues.
              global_nuclei.push_back(nuc);
              ConditionalOStreams::pout_base()
                << "  New nucleus at: " << nuc.location << "\n";
              ++count;
              any_nucleation_occurred = true;
            }
          new_nuclei.pop_back();
        }
      ConditionalOStreams::pout_base() << "  " << count
                                       << " new nuclei after exclusion.\n"
                                          "  "
                                       << global_nuclei.size() << " total nuclei.\n\n"
                                       << std::flush;
    }
  MPI_Bcast(&any_nucleation_occurred, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
  mpi_broadcast_nuclei(global_nuclei);
  return any_nucleation_occurred;
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
NucleationManager<dim, degree, number>::mpi_gather_nuclei(
  std::vector<Nucleus<dim>> &local_nuclei)
{
  // Step 1: Share how many nuclei each rank has
  int              rank        = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  int              num_procs   = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  int              local_count = local_nuclei.size();
  std::vector<int> nuclei_counts_per_rank;
  if (rank == 0)
    {
      nuclei_counts_per_rank.resize(num_procs);
    }

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
  if (rank == 0)
    {
      local_nuclei = gathered_nuclei;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
NucleationManager<dim, degree, number>::mpi_broadcast_nuclei(
  std::vector<Nucleus<dim>> &local_nuclei)
{
  int rank  = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  int count = local_nuclei.size();
  MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0)
    {
      local_nuclei.resize(count);
    }

  MPI_Bcast(local_nuclei.data(), count, Nucleus<dim>::mpi_datatype(), 0, MPI_COMM_WORLD);
}

PRISMS_PF_END_NAMESPACE
