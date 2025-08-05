// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/lac/solver_cg.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/linear_solver_base.h>
#include <prismspf/solvers/linear_solver_identity.h>

#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
IdentitySolver<dim, degree, number>::IdentitySolver(
  const SolverContext<dim, degree, number> &_solver_context,
  const VariableAttributes                 &_variable_attributes)
  : LinearSolverBase<dim, degree, number>(_solver_context, _variable_attributes)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
IdentitySolver<dim, degree, number>::init()
{
  // Call the base class init
  this->LinearSolverBase<dim, degree, number>::init();

  // Add some stuff to the matrix free operator
  this->get_system_matrix()->clear();
  this->get_system_matrix()->initialize(
    this->get_matrix_free_container().get_matrix_free(),
    this->get_element_volume_container().get_element_volume());
  this->get_update_system_matrix()->clear();
  this->get_update_system_matrix()->initialize(
    this->get_matrix_free_container().get_matrix_free(),
    this->get_element_volume_container().get_element_volume());

  this->get_system_matrix()->add_global_to_local_mapping(
    this->get_residual_global_to_local_solution());
  this->get_system_matrix()->add_src_solution_subset(this->get_residual_src());

  this->get_update_system_matrix()->add_global_to_local_mapping(
    this->get_newton_update_global_to_local_solution());
  this->get_update_system_matrix()->add_src_solution_subset(
    this->get_newton_update_src());

  // Apply constraints
  this->get_constraint_handler()
    .get_constraint(this->get_field_index())
    .distribute(
      *(this->get_solution_handler().get_solution_vector(this->get_field_index(),
                                                         DependencyType::Normal)));
}

template <unsigned int dim, unsigned int degree, typename number>
void
IdentitySolver<dim, degree, number>::reinit()
{
  // Call the base class reinit
  this->LinearSolverBase<dim, degree, number>::reinit();

  // Add some stuff to the matrix free operator
  this->get_system_matrix()->clear();
  this->get_system_matrix()->initialize(
    this->get_matrix_free_container().get_matrix_free(),
    this->get_element_volume_container().get_element_volume());
  this->get_update_system_matrix()->clear();
  this->get_update_system_matrix()->initialize(
    this->get_matrix_free_container().get_matrix_free(),
    this->get_element_volume_container().get_element_volume());

  this->get_system_matrix()->add_global_to_local_mapping(
    this->get_residual_global_to_local_solution());
  this->get_system_matrix()->add_src_solution_subset(this->get_residual_src());

  this->get_update_system_matrix()->add_global_to_local_mapping(
    this->get_newton_update_global_to_local_solution());
  this->get_update_system_matrix()->add_src_solution_subset(
    this->get_newton_update_src());

  // Apply constraints
  this->get_constraint_handler()
    .get_constraint(this->get_field_index())
    .distribute(
      *(this->get_solution_handler().get_solution_vector(this->get_field_index(),
                                                         DependencyType::Normal)));
}

template <unsigned int dim, unsigned int degree, typename number>
void
IdentitySolver<dim, degree, number>::solve(const number &step_length)
{
  auto *solution =
    this->get_solution_handler().get_solution_vector(this->get_field_index(),
                                                     DependencyType::Normal);

  // Compute the residual
  this->get_system_matrix()->compute_residual(*this->get_residual(), *solution);
  if (this->get_user_inputs().get_output_parameters().should_output(
        this->get_user_inputs().get_temporal_discretization().get_increment()))
    {
      ConditionalOStreams::pout_summary()
        << "  field: " << this->get_field_index()
        << " Initial residual: " << this->get_residual()->l2_norm() << std::flush;
    }

  ConditionalOStreams::pout_summary()
    << "\n Field: " << this->get_field_index()
    << " Initial Residual: " << this->get_residual()->l2_norm() << " size "
    << this->get_residual()->size()
    << " Initial Newton Update: " << this->get_newton_update()->l2_norm() << " size "
    << this->get_newton_update()->size() << std::flush;

  // Determine the residual tolerance
  this->compute_solver_tolerance();

  // Update solver controls
  this->get_solver_control().set_tolerance(this->get_tolerance());
  dealii::SolverCG<VectorType> cg_solver(this->get_solver_control());

  try
    {
      *this->get_newton_update() = 0.0;
      cg_solver.solve(*(this->get_update_system_matrix()),
                      *(this->get_newton_update()),
                      *(this->get_residual()),
                      dealii::PreconditionIdentity());
    }
  catch (...)
    {
      ConditionalOStreams::pout_base()
        << "Warning: linear solver did not converge as per set tolerances.\n";
    }
  this->get_constraint_handler()
    .get_constraint(this->get_field_index())
    .set_zero(*this->get_newton_update());

  if (this->get_user_inputs().get_output_parameters().should_output(
        this->get_user_inputs().get_temporal_discretization().get_increment()))
    {
      ConditionalOStreams::pout_summary()
        << " Final residual: " << this->get_solver_control().last_value()
        << " Steps: " << this->get_solver_control().last_step() << "\n"
        << std::flush;
    }

  // Update the solutions
  (*solution).add(step_length, *this->get_newton_update());

  // Apply constraints
  // This may be redundant with the constraints on the update step.
  this->get_constraint_handler()
    .get_constraint(this->get_field_index())
    .distribute(*solution);
}

#include "solvers/linear_solver_identity.inst"

PRISMS_PF_END_NAMESPACE
