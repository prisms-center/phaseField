// solveIncrement() method for MatrixFreePDE class

#include <deal.II/lac/solver_cg.h>

#include "../../include/matrixFreePDE.h"

// solve each time increment
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::solveIncrement(bool skip_time_dependent)
{
  // log time
  computing_timer.enter_subsection("matrixFreePDE: solveIncrements");
  Timer time;
  char  buffer[200];

  // Get the RHS of the explicit equations
  if (hasExplicitEquation && !skip_time_dependent)
    {
      computeExplicitRHS();
    }

  // solve for each field
  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      currentFieldIndex = fieldIndex; // Used in computeLHS()

      // Parabolic (first order derivatives in time) fields
      if (fields[fieldIndex].pdetype == EXPLICIT_TIME_DEPENDENT && !skip_time_dependent)
        {
          updateExplicitSolution(fieldIndex);

          // Apply Boundary conditions
          applyBCs(fieldIndex);

          // Print update to screen and confirm that solution isn't nan
          if (currentIncrement % userInputs.skip_print_steps == 0)
            {
              double solution_L2_norm = solutionSet[fieldIndex]->l2_norm();

              snprintf(buffer,
                       sizeof(buffer),
                       "field '%2s' [explicit solve]: current solution: "
                       "%12.6e, current residual:%12.6e\n",
                       fields[fieldIndex].name.c_str(),
                       solution_L2_norm,
                       residualSet[fieldIndex]->l2_norm());
              pcout << buffer;

              if (!numbers::is_finite(solution_L2_norm))
                {
                  snprintf(buffer,
                           sizeof(buffer),
                           "ERROR: field '%s' solution is NAN. exiting.\n\n",
                           fields[fieldIndex].name.c_str());
                  pcout << buffer;
                  exit(-1);
                }
            }
        }
    }

  // Now, update the non-explicit variables
  // For the time being, this is just the elliptic equations, but implicit
  // parabolic and auxilary equations should also be here
  if (hasNonExplicitEquation)
    {
      bool         nonlinear_it_converged = false;
      unsigned int nonlinear_it_index     = 0;

      while (!nonlinear_it_converged)
        {
          nonlinear_it_converged = true; // Set to true here and will be set to false if
                                         // any variable isn't converged

          // Update residualSet for the non-explicitly updated variables
          // compute_nonexplicit_RHS()
          // Ideally, I'd just do this for the non-explicit variables, but for
          // now I'll do all of them this is a little redundant, but hopefully
          // not too terrible
          computeNonexplicitRHS();

          for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
            {
              currentFieldIndex = fieldIndex; // Used in computeLHS()

              if ((fields[fieldIndex].pdetype == IMPLICIT_TIME_DEPENDENT &&
                   !skip_time_dependent) ||
                  fields[fieldIndex].pdetype == TIME_INDEPENDENT)
                {
                  if (currentIncrement % userInputs.skip_print_steps == 0 &&
                      userInputs.var_nonlinear[fieldIndex])
                    {
                      snprintf(buffer,
                               sizeof(buffer),
                               "field '%2s' [nonlinear solve]: current "
                               "solution: %12.6e, current residual:%12.6e\n",
                               fields[fieldIndex].name.c_str(),
                               solutionSet[fieldIndex]->l2_norm(),
                               residualSet[fieldIndex]->l2_norm());
                      pcout << buffer;
                    }

                  nonlinear_it_converged =
                    updateImplicitSolution(fieldIndex, nonlinear_it_index);

                  // Apply Boundary conditions
                  applyBCs(fieldIndex);
                }
              else if (fields[fieldIndex].pdetype == AUXILIARY)
                {
                  if (userInputs.var_nonlinear[fieldIndex] || nonlinear_it_index == 0)
                    {
                      // If the equation for this field is nonlinear, save the old
                      // solution
                      if (userInputs.var_nonlinear[fieldIndex])
                        {
                          if (fields[fieldIndex].type == SCALAR)
                            {
                              dU_scalar = *solutionSet[fieldIndex];
                            }
                          else
                            {
                              dU_vector = *solutionSet[fieldIndex];
                            }
                        }

                      updateExplicitSolution(fieldIndex);

                      // Apply Boundary conditions
                      applyBCs(fieldIndex);

                      // Print update to screen
                      if (currentIncrement % userInputs.skip_print_steps == 0)
                        {
                          snprintf(buffer,
                                   sizeof(buffer),
                                   "field '%2s' [auxiliary solve]: current solution: "
                                   "%12.6e, current residual:%12.6e\n",
                                   fields[fieldIndex].name.c_str(),
                                   solutionSet[fieldIndex]->l2_norm(),
                                   residualSet[fieldIndex]->l2_norm());
                          pcout << buffer;
                        }

                      // Check to see if this individual variable has converged
                      if (userInputs.var_nonlinear[fieldIndex])
                        {
                          if (userInputs.nonlinear_solver_parameters.getToleranceType(
                                fieldIndex) == ABSOLUTE_SOLUTION_CHANGE)
                            {
                              double diff;

                              if (fields[fieldIndex].type == SCALAR)
                                {
                                  dU_scalar -= *solutionSet[fieldIndex];
                                  diff = dU_scalar.l2_norm();
                                }
                              else
                                {
                                  dU_vector -= *solutionSet[fieldIndex];
                                  diff = dU_vector.l2_norm();
                                }
                              if (currentIncrement % userInputs.skip_print_steps == 0)
                                {
                                  pcout << "Relative difference between nonlinear "
                                           "iterations: "
                                        << diff << " " << nonlinear_it_index << " "
                                        << currentIncrement << std::endl;
                                }

                              if (diff > userInputs.nonlinear_solver_parameters
                                           .getToleranceValue(fieldIndex) &&
                                  nonlinear_it_index <
                                    userInputs.nonlinear_solver_parameters
                                      .getMaxIterations())
                                {
                                  nonlinear_it_converged = false;
                                }
                            }
                          else
                            {
                              std::cerr << "PRISMS-PF Error: Nonlinear solver tolerance "
                                           "types other than ABSOLUTE_CHANGE have yet to "
                                           "be implemented."
                                        << std::endl;
                            }
                        }
                    }
                }

              // check if solution is nan
              if (!numbers::is_finite(solutionSet[fieldIndex]->l2_norm()))
                {
                  snprintf(buffer,
                           sizeof(buffer),
                           "ERROR: field '%s' solution is NAN. exiting.\n\n",
                           fields[fieldIndex].name.c_str());
                  pcout << buffer;
                  exit(-1);
                }
            }

          nonlinear_it_index++;
        }
    }

  if (currentIncrement % userInputs.skip_print_steps == 0)
    {
      pcout << "wall time: " << time.wall_time() << "s\n";
    }
  // log time
  computing_timer.leave_subsection("matrixFreePDE: solveIncrements");
}

// Application of boundary conditions
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::applyBCs(unsigned int fieldIndex)
{
  // Add Neumann BCs
  if (fields[fieldIndex].hasNeumannBCs)
    {
      // Currently commented out because it isn't working yet
      // applyNeumannBCs();
    }

  // Set the Dirichelet values (hanging node constraints don't need to be distributed
  // every time step, only at output)
  if (fields[fieldIndex].hasDirichletBCs)
    {
      // Apply non-uniform Dirlichlet_BCs to the current field
      if (fields[fieldIndex].hasnonuniformDirichletBCs)
        {
          DoFHandler<dim> *dof_handler;
          dof_handler = dofHandlersSet_nonconst.at(currentFieldIndex);
          IndexSet *locally_relevant_dofs;
          locally_relevant_dofs = locally_relevant_dofsSet_nonconst.at(currentFieldIndex);
          locally_relevant_dofs->clear();
          DoFTools::extract_locally_relevant_dofs(*dof_handler, *locally_relevant_dofs);
          AffineConstraints<double> *constraintsDirichlet;
          constraintsDirichlet = constraintsDirichletSet_nonconst.at(currentFieldIndex);
          constraintsDirichlet->clear();
          constraintsDirichlet->reinit(*locally_relevant_dofs);
          applyDirichletBCs();
          constraintsDirichlet->close();
        }
      // Distribute for Uniform or Non-Uniform Dirichlet BCs
      constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
    }
  solutionSet[fieldIndex]->update_ghost_values();
}

// Explicit time step for matrixfree solve
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::updateExplicitSolution(unsigned int fieldIndex)
{
  // Explicit-time step each DOF
  // Takes advantage of knowledge that the length of solutionSet and residualSet
  // is an integer multiple of the length of invM for vector variables
  if (fields[fieldIndex].type == SCALAR)
    {
#if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR < 4)
      unsigned int invM_size = invMscalar.local_size();
      for (unsigned int dof = 0; dof < solutionSet[fieldIndex]->local_size(); ++dof)
        {
#else
      unsigned int invM_size = invMscalar.locally_owned_size();
      for (unsigned int dof = 0; dof < solutionSet[fieldIndex]->locally_owned_size();
           ++dof)
        {
#endif
          solutionSet[fieldIndex]->local_element(dof) =
            invMscalar.local_element(dof % invM_size) *
            residualSet[fieldIndex]->local_element(dof);
        }
    }
  else if (fields[fieldIndex].type == VECTOR)
    {
#if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR < 4)
      unsigned int invM_size = invMvector.local_size();
      for (unsigned int dof = 0; dof < solutionSet[fieldIndex]->local_size(); ++dof)
        {
#else
      unsigned int invM_size = invMvector.locally_owned_size();
      for (unsigned int dof = 0; dof < solutionSet[fieldIndex]->locally_owned_size();
           ++dof)
        {
#endif
          solutionSet[fieldIndex]->local_element(dof) =
            invMvector.local_element(dof % invM_size) *
            residualSet[fieldIndex]->local_element(dof);
        }
    }
}

template <int dim, int degree>
bool
MatrixFreePDE<dim, degree>::updateImplicitSolution(unsigned int fieldIndex,
                                                   unsigned int nonlinear_it_index)
{
  char buffer[200];

  // Assume convergence criterion is met, unless otherwise proven later on.
  bool nonlinear_it_converged = true;

  // Apply Dirichlet BC's. This clears the residual where we want to apply Dirichlet BCs,
  // otherwise the solver sees a positive residual
  constraintsDirichletSet[fieldIndex]->set_zero(*residualSet[fieldIndex]);

  // Grab solver controls
  double tol_value;
  if (userInputs.linear_solver_parameters.getToleranceType(fieldIndex) ==
      ABSOLUTE_RESIDUAL)
    {
      tol_value = userInputs.linear_solver_parameters.getToleranceValue(fieldIndex);
    }
  else
    {
      tol_value = userInputs.linear_solver_parameters.getToleranceValue(fieldIndex) *
                  residualSet[fieldIndex]->l2_norm();
    }

  SolverControl solver_control(userInputs.linear_solver_parameters.getMaxIterations(
                                 fieldIndex),
                               tol_value);

  // Currently the only allowed solver is SolverCG, the
  // SolverType input variable is a dummy
  SolverCG<vectorType> solver(solver_control);

  // Solve
  try
    {
      if (fields[fieldIndex].type == SCALAR)
        {
          dU_scalar = 0.0;
          solver.solve(*this,
                       dU_scalar,
                       *residualSet[fieldIndex],
                       IdentityMatrix(solutionSet[fieldIndex]->size()));
        }
      else
        {
          dU_vector = 0.0;
          solver.solve(*this,
                       dU_vector,
                       *residualSet[fieldIndex],
                       IdentityMatrix(solutionSet[fieldIndex]->size()));
        }
    }
  catch (...)
    {
      pcout << "\nWarning: linear solver did not converge as "
               "per set tolerances. consider increasing the "
               "maximum number of iterations or decreasing the "
               "solver tolerance.\n";
    }

  if (userInputs.var_nonlinear[fieldIndex])
    {
      // Now that we have the calculated change in the solution,
      // we need to select a damping coefficient
      double damping_coefficient;

      if (userInputs.nonlinear_solver_parameters.getBacktrackDampingFlag(fieldIndex))
        {
          vectorType solutionSet_old = *solutionSet[fieldIndex];
          double     residual_old    = residualSet[fieldIndex]->l2_norm();

          damping_coefficient            = 1.0;
          bool damping_coefficient_found = false;
          while (!damping_coefficient_found)
            {
              if (fields[fieldIndex].type == SCALAR)
                {
                  solutionSet[fieldIndex]->sadd(1.0, damping_coefficient, dU_scalar);
                }
              else
                {
                  solutionSet[fieldIndex]->sadd(1.0, damping_coefficient, dU_vector);
                }

              computeNonexplicitRHS();

              for (const auto &it : *valuesDirichletSet[fieldIndex])
                {
                  if (residualSet[fieldIndex]->in_local_range(it.first))
                    {
                      (*residualSet[fieldIndex])(it.first) = 0.0;
                    }
                }

              double residual_new = residualSet[fieldIndex]->l2_norm();

              if (currentIncrement % userInputs.skip_print_steps == 0)
                {
                  pcout << "    Old residual: " << residual_old
                        << " Damping Coeff: " << damping_coefficient
                        << " New Residual: " << residual_new << std::endl;
                }

              // An improved approach would use the
              // Armijo–Goldstein condition to ensure a
              // sufficent decrease in the residual. This way is
              // just scales the residual.
              if ((residual_new <
                   (residual_old * userInputs.nonlinear_solver_parameters
                                     .getBacktrackResidualDecreaseCoeff(fieldIndex))) ||
                  damping_coefficient < 1.0e-4)
                {
                  damping_coefficient_found = true;
                }
              else
                {
                  damping_coefficient *=
                    userInputs.nonlinear_solver_parameters.getBacktrackStepModifier(
                      fieldIndex);
                  *solutionSet[fieldIndex] = solutionSet_old;
                }
            }
        }
      else
        {
          damping_coefficient =
            userInputs.nonlinear_solver_parameters.getDefaultDampingCoefficient(
              fieldIndex);

          if (fields[fieldIndex].type == SCALAR)
            {
              solutionSet[fieldIndex]->sadd(1.0, damping_coefficient, dU_scalar);
            }
          else
            {
              solutionSet[fieldIndex]->sadd(1.0, damping_coefficient, dU_vector);
            }
        }

      if (currentIncrement % userInputs.skip_print_steps == 0)
        {
          double dU_norm;
          if (fields[fieldIndex].type == SCALAR)
            {
              dU_norm = dU_scalar.l2_norm();
            }
          else
            {
              dU_norm = dU_vector.l2_norm();
            }
          snprintf(buffer,
                   sizeof(buffer),
                   "field '%2s' [linear solve]: initial "
                   "residual:%12.6e, current residual:%12.6e, "
                   "nsteps:%u, tolerance criterion:%12.6e, "
                   "solution: %12.6e, dU: %12.6e\n",
                   fields[fieldIndex].name.c_str(),
                   residualSet[fieldIndex]->l2_norm(),
                   solver_control.last_value(),
                   solver_control.last_step(),
                   solver_control.tolerance(),
                   solutionSet[fieldIndex]->l2_norm(),
                   dU_norm);
          pcout << buffer;
        }

      // Check to see if this individual variable has converged
      if (userInputs.nonlinear_solver_parameters.getToleranceType(fieldIndex) ==
          ABSOLUTE_SOLUTION_CHANGE)
        {
          double diff;

          if (fields[fieldIndex].type == SCALAR)
            {
              diff = dU_scalar.l2_norm();
            }
          else
            {
              diff = dU_vector.l2_norm();
            }
          if (currentIncrement % userInputs.skip_print_steps == 0)
            {
              pcout << "Relative difference between nonlinear "
                       "iterations: "
                    << diff << " " << nonlinear_it_index << " " << currentIncrement
                    << std::endl;
            }

          if (diff >
                userInputs.nonlinear_solver_parameters.getToleranceValue(fieldIndex) &&
              nonlinear_it_index <
                userInputs.nonlinear_solver_parameters.getMaxIterations())
            {
              nonlinear_it_converged = false;
            }
        }
      else
        {
          std::cerr << "PRISMS-PF Error: Nonlinear solver tolerance "
                       "types other than ABSOLUTE_CHANGE have yet to "
                       "be implemented."
                    << std::endl;
        }
    }
  else
    {
      if (nonlinear_it_index == 0)
        {
          if (fields[fieldIndex].type == SCALAR)
            {
              *solutionSet[fieldIndex] += dU_scalar;
            }
          else
            {
              *solutionSet[fieldIndex] += dU_vector;
            }

          if (currentIncrement % userInputs.skip_print_steps == 0)
            {
              double dU_norm;
              if (fields[fieldIndex].type == SCALAR)
                {
                  dU_norm = dU_scalar.l2_norm();
                }
              else
                {
                  dU_norm = dU_vector.l2_norm();
                }
              snprintf(buffer,
                       sizeof(buffer),
                       "field '%2s' [linear solve]: initial "
                       "residual:%12.6e, current residual:%12.6e, "
                       "nsteps:%u, tolerance criterion:%12.6e, "
                       "solution: %12.6e, dU: %12.6e\n",
                       fields[fieldIndex].name.c_str(),
                       residualSet[fieldIndex]->l2_norm(),
                       solver_control.last_value(),
                       solver_control.last_step(),
                       solver_control.tolerance(),
                       solutionSet[fieldIndex]->l2_norm(),
                       dU_norm);
              pcout << buffer;
            }
        }
    }

  return nonlinear_it_converged;
}

#include "../../include/matrixFreePDE_template_instantiations.h"
