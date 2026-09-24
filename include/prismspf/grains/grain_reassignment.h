// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mpi.h>
#include <deal.II/fe/fe_values.h>

#include <prismspf/core/pde_operator_base.h>
#include <prismspf/core/simulation_timer.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/types.h>

#include <prismspf/grains/flood_filler.h>
#include <prismspf/grains/grains.h>

#include <prismspf/solvers/solve_context.h>

#include <prismspf/user_inputs/miscellaneous_parameters.h>
#include <prismspf/user_inputs/temporal_discretization.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/timer.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

// TODO: A lot of stuff in here could be static and in namespaces instead of classes.

/**
 * This class contains methods to make changes to lists of
 * SimplifiedGrainRepresentation objects to aid in order parameter remapping.
 * Currently, this is a stub class containing two unrelated methods and no
 * internal members. If this stays the case, this class may be absorbed into
 * another one.
 */
template <unsigned int dim>
class SimplifiedGrainManipulator
{
public:
  /**
   * This method checks for collisions between SimplifiedGrainRepresentation
   * objects with the same order parameter and reassigns them, if needed.
   */
  static void
  reassign_grains(std::vector<SimplifiedGrainRepresentation<dim>> &grain_representations,
                  double                                           buffer_distance,
                  int                                 number_of_remapped_fields,
                  const std::vector<FieldAttributes> &field_attributes);

  /**
   * This method checks the centers of two lists of
   * SimplifiedGrainRepresentation objects from different times in the
   * simulation and reassigns the grain ids so that they consistently refer to
   * the same grains.
   */
  static void
  transfer_grain_ids(
    const std::vector<SimplifiedGrainRepresentation<dim>> &old_grain_representations,
    std::vector<SimplifiedGrainRepresentation<dim>>       &new_grain_representations);
};

/**
 * This class uses information from the list of SimplifiedGrainRepresentation
 * objects to reassign grains across multiple solution fields.
 */
template <unsigned int dim, typename number>
class OrderParameterRemapper
{
public:
  /**
   * This method does the core work of the class to reassign grains across
   * solution vectors based on the list of SimplifiedGrainRepresentation
   * objects.
   */
  static void
  remap(std::vector<SimplifiedGrainRepresentation<dim>> &grain_representations,
        SolutionIndexer<dim, number>                    &solution_fields,
        const dealii::DoFHandler<dim>                   &dof_handler,
        unsigned int                                     dofs_per_cell);
};

template <unsigned int dim, unsigned int degree, typename number>
class GrainReassignmentManager
{
public:
  static void
  reassign_grains(
    const SolveContext<dim, degree, number>         &solve_context,
    std::vector<SimplifiedGrainRepresentation<dim>> &simplified_grain_representations);
};

PRISMS_PF_END_NAMESPACE
