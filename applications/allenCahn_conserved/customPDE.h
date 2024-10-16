#include "matrixFreePDE.h"

using namespace dealii;

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
  // Constructor
  customPDE(userInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
    , userInputs(_userInputs) {};

  // Function to set the initial conditions (in ICs_and_BCs.h)
  void
  setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                      [[maybe_unused]] const unsigned int index,
                      [[maybe_unused]] double            &scalar_IC,
                      [[maybe_unused]] Vector<double>    &vector_IC) override;

  // Function to set the non-uniform Dirichlet boundary conditions (in
  // ICs_and_BCs.h)
  void
  setNonUniformDirichletBCs([[maybe_unused]] const Point<dim>  &p,
                            [[maybe_unused]] const unsigned int index,
                            [[maybe_unused]] const unsigned int direction,
                            [[maybe_unused]] const double       time,
                            [[maybe_unused]] double            &scalar_BC,
                            [[maybe_unused]] Vector<double>    &vector_BC) override;

private:
#include "typeDefs.h"

  const userInputParameters<dim> userInputs;

  // Function to set the RHS of the governing equations for explicit time
  // dependent equations (in equations.cc)
  void
  explicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                        &variable_list,
    [[maybe_unused]] Point<dim, VectorizedArray<double>> q_point_loc) const override;

  // Function to set the RHS of the governing equations for all other equations
  // (in equations.h)
  void
  nonExplicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                        &variable_list,
    [[maybe_unused]] Point<dim, VectorizedArray<double>> q_point_loc) const override;

  // Function to set the LHS of the governing equations (in equations.cc)
  void
  equationLHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                        &variable_list,
    [[maybe_unused]] Point<dim, VectorizedArray<double>> q_point_loc) const override;

// Function to set postprocessing expressions (in postprocess.cc)
#ifdef POSTPROCESS_FILE_EXISTS
  void
  postProcessedFields(
    [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
      &variable_list,
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &pp_variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc)
    const override;
#endif

// Function to set the nucleation probability (in nucleation.cc)
#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability([[maybe_unused]] variableValueContainer variable_value,
                           [[maybe_unused]] double                 dV) const override;
#endif

  // ================================================================
  // Methods specific to this subclass
  // ================================================================

  // Function to override solveIncrement from
  // ../../src/matrixfree/solveIncrement.cc
  void
  solveIncrement(bool skip_time_dependent) override;

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  double MnV = userInputs.get_model_constant_double("MnV");
  double KnV = userInputs.get_model_constant_double("KnV");

  double integrated_mu;
  double integrated_n;

  // ================================================================
};

// =================================================================================
// Function overriding solveIncrement
// =================================================================================
// solve each time increment
#include <deal.II/lac/solver_cg.h>

template <int dim, int degree>
void
customPDE<dim, degree>::solveIncrement(bool skip_time_dependent)
{
  // log time
  this->computing_timer.enter_subsection("matrixFreePDE: solveIncrements");
  Timer time;
  char  buffer[200];

  // Calculating integral for mu (field 1)
  this->computeIntegral(integrated_n, 0, this->solutionSet);
  this->computeIntegral(integrated_mu, 1, this->solutionSet);

  if (this->currentIncrement % userInputs.skip_print_steps == 0)
    {
      snprintf(buffer, sizeof(buffer), "Integrated mu is %12.6e\n", integrated_mu);
      this->pcout << buffer;
    }

  // Get the RHS of the explicit equations
  if (this->hasExplicitEquation && !skip_time_dependent)
    {
      this->computeExplicitRHS();
    }

  // solve for each field
  for (unsigned int fieldIndex = 0; fieldIndex < this->fields.size(); fieldIndex++)
    {
      this->currentFieldIndex = fieldIndex; // Used in computeLHS()

      // Parabolic (first order derivatives in time) fields
      if (this->fields[fieldIndex].pdetype == EXPLICIT_TIME_DEPENDENT &&
          !skip_time_dependent)
        {
          this->updateExplicitSolution(fieldIndex);

          // Apply Boundary conditions
          this->applyBCs(fieldIndex);

          // Print update to screen and confirm that solution isn't nan
          if (this->currentIncrement % userInputs.skip_print_steps == 0)
            {
              double solution_L2_norm = this->solutionSet[fieldIndex]->l2_norm();

              snprintf(buffer,
                       sizeof(buffer),
                       "field '%2s' [explicit solve]: current solution: "
                       "%12.6e, current residual:%12.6e\n",
                       this->fields[fieldIndex].name.c_str(),
                       solution_L2_norm,
                       this->residualSet[fieldIndex]->l2_norm());
              this->pcout << buffer;

              if (!numbers::is_finite(solution_L2_norm))
                {
                  snprintf(buffer,
                           sizeof(buffer),
                           "ERROR: field '%s' solution is NAN. exiting.\n\n",
                           this->fields[fieldIndex].name.c_str());
                  this->pcout << buffer;
                  exit(-1);
                }
            }
        }
    }

  // Now, update the non-explicit variables
  // For the time being, this is just the elliptic equations, but implicit
  // parabolic and auxilary equations should also be here
  if (this->hasNonExplicitEquation)
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
          this->computeNonexplicitRHS();

          for (unsigned int fieldIndex = 0; fieldIndex < this->fields.size();
               fieldIndex++)
            {
              this->currentFieldIndex = fieldIndex; // Used in computeLHS()

              if ((this->fields[fieldIndex].pdetype == IMPLICIT_TIME_DEPENDENT &&
                   !skip_time_dependent) ||
                  this->fields[fieldIndex].pdetype == TIME_INDEPENDENT)
                {
                  if (this->currentIncrement % userInputs.skip_print_steps == 0 &&
                      userInputs.var_nonlinear[fieldIndex])
                    {
                      snprintf(buffer,
                               sizeof(buffer),
                               "field '%2s' [nonlinear solve]: current "
                               "solution: %12.6e, current residual:%12.6e\n",
                               this->fields[fieldIndex].name.c_str(),
                               this->solutionSet[fieldIndex]->l2_norm(),
                               this->residualSet[fieldIndex]->l2_norm());
                      this->pcout << buffer;
                    }

                  LinearAlgebra::distributed::Vector<double> solution_diff =
                    *this->solutionSet[fieldIndex];

                  // apply Dirichlet BC's
                  //  This clears the residual where we want to apply Dirichlet
                  //  BCs, otherwise the solver sees a positive residual
                  this->constraintsDirichletSet[fieldIndex]->set_zero(
                    *this->residualSet[fieldIndex]);

                  // solver controls
                  double tol_value;
                  if (MatrixFreePDE<dim, degree>::userInputs.linear_solver_parameters
                        .getToleranceType(fieldIndex) == ABSOLUTE_RESIDUAL)
                    {
                      tol_value =
                        MatrixFreePDE<dim, degree>::userInputs.linear_solver_parameters
                          .getToleranceValue(fieldIndex);
                    }
                  else
                    {
                      tol_value =
                        MatrixFreePDE<dim, degree>::userInputs.linear_solver_parameters
                          .getToleranceValue(fieldIndex) *
                        this->residualSet[fieldIndex]->l2_norm();
                    }

                  SolverControl solver_control(
                    MatrixFreePDE<dim, degree>::userInputs.linear_solver_parameters
                      .getMaxIterations(fieldIndex),
                    tol_value);

                  // Currently the only allowed solver is SolverCG, the
                  // SolverType input variable is a dummy
                  SolverCG<vectorType> solver(solver_control);

                  // solve
                  try
                    {
                      if (this->fields[fieldIndex].type == SCALAR)
                        {
                          this->dU_scalar = 0.0;
                          solver.solve(*this,
                                       this->dU_scalar,
                                       *this->residualSet[fieldIndex],
                                       IdentityMatrix(
                                         this->solutionSet[fieldIndex]->size()));
                        }
                      else
                        {
                          this->dU_vector = 0.0;
                          solver.solve(*this,
                                       this->dU_vector,
                                       *this->residualSet[fieldIndex],
                                       IdentityMatrix(
                                         this->solutionSet[fieldIndex]->size()));
                        }
                    }
                  catch (...)
                    {
                      this->pcout << "\nWarning: linear solver did not converge as per "
                                     "set tolerances. consider increasing the maximum "
                                     "number of iterations or decreasing the solver "
                                     "tolerance.\n";
                    }

                  if (userInputs.var_nonlinear[fieldIndex])
                    {
                      // Now that we have the calculated change in the solution,
                      // we need to select a damping coefficient
                      double damping_coefficient;

                      if (MatrixFreePDE<dim, degree>::userInputs
                            .nonlinear_solver_parameters.getBacktrackDampingFlag(
                              fieldIndex))
                        {
                          vectorType solutionSet_old = *this->solutionSet[fieldIndex];
                          double residual_old = this->residualSet[fieldIndex]->l2_norm();

                          damping_coefficient            = 1.0;
                          bool damping_coefficient_found = false;
                          while (!damping_coefficient_found)
                            {
                              if (this->fields[fieldIndex].type == SCALAR)
                                {
                                  this->solutionSet[fieldIndex]->sadd(1.0,
                                                                      damping_coefficient,
                                                                      this->dU_scalar);
                                }
                              else
                                {
                                  this->solutionSet[fieldIndex]->sadd(1.0,
                                                                      damping_coefficient,
                                                                      this->dU_vector);
                                }

                              this->computeNonexplicitRHS();

                              for (const auto &it : *this->valuesDirichletSet[fieldIndex])
                                {
                                  if (this->residualSet[fieldIndex]->in_local_range(
                                        it.first))
                                    {
                                      (*this->residualSet[fieldIndex])(it.first) = 0.0;
                                    }
                                }

                              double residual_new =
                                this->residualSet[fieldIndex]->l2_norm();

                              if (this->currentIncrement % userInputs.skip_print_steps ==
                                  0)
                                {
                                  this->pcout << "    Old residual: " << residual_old
                                              << " Damping Coeff: " << damping_coefficient
                                              << " New Residual: " << residual_new
                                              << std::endl;
                                }

                              // An improved approach would use the
                              // Armijo–Goldstein condition to ensure a
                              // sufficent decrease in the residual. This way is
                              // just scales the residual.
                              if ((residual_new <
                                   (residual_old *
                                    MatrixFreePDE<dim, degree>::userInputs
                                      .nonlinear_solver_parameters
                                      .getBacktrackResidualDecreaseCoeff(fieldIndex))) ||
                                  damping_coefficient < 1.0e-4)
                                {
                                  damping_coefficient_found = true;
                                }
                              else
                                {
                                  damping_coefficient *=
                                    MatrixFreePDE<dim, degree>::userInputs
                                      .nonlinear_solver_parameters
                                      .getBacktrackStepModifier(fieldIndex);
                                  *this->solutionSet[fieldIndex] = solutionSet_old;
                                }
                            }
                        }
                      else
                        {
                          damping_coefficient =
                            MatrixFreePDE<dim, degree>::userInputs
                              .nonlinear_solver_parameters.getDefaultDampingCoefficient(
                                fieldIndex);

                          if (this->fields[fieldIndex].type == SCALAR)
                            {
                              this->solutionSet[fieldIndex]->sadd(1.0,
                                                                  damping_coefficient,
                                                                  this->dU_scalar);
                            }
                          else
                            {
                              this->solutionSet[fieldIndex]->sadd(1.0,
                                                                  damping_coefficient,
                                                                  this->dU_vector);
                            }
                        }

                      if (this->currentIncrement % userInputs.skip_print_steps == 0)
                        {
                          double dU_norm;
                          if (this->fields[fieldIndex].type == SCALAR)
                            {
                              dU_norm = this->dU_scalar.l2_norm();
                            }
                          else
                            {
                              dU_norm = this->dU_vector.l2_norm();
                            }
                          snprintf(buffer,
                                   sizeof(buffer),
                                   "field '%2s' [linear solve]: initial "
                                   "residual:%12.6e, current residual:%12.6e, "
                                   "nsteps:%u, tolerance criterion:%12.6e, "
                                   "solution: %12.6e, dU: %12.6e\n",
                                   this->fields[fieldIndex].name.c_str(),
                                   this->residualSet[fieldIndex]->l2_norm(),
                                   solver_control.last_value(),
                                   solver_control.last_step(),
                                   solver_control.tolerance(),
                                   this->solutionSet[fieldIndex]->l2_norm(),
                                   dU_norm);
                          this->pcout << buffer;
                        }

                      // Check to see if this individual variable has converged
                      if (MatrixFreePDE<dim, degree>::userInputs
                            .nonlinear_solver_parameters.getToleranceType(fieldIndex) ==
                          ABSOLUTE_SOLUTION_CHANGE)
                        {
                          double diff;

                          if (this->fields[fieldIndex].type == SCALAR)
                            {
                              diff = this->dU_scalar.l2_norm();
                            }
                          else
                            {
                              diff = this->dU_vector.l2_norm();
                            }
                          if (this->currentIncrement % userInputs.skip_print_steps == 0)
                            {
                              this->pcout << "Relative difference between "
                                             "nonlinear iterations: "
                                          << diff << " " << nonlinear_it_index << " "
                                          << this->currentIncrement << std::endl;
                            }

                          if (diff > MatrixFreePDE<dim, degree>::userInputs
                                       .nonlinear_solver_parameters.getToleranceValue(
                                         fieldIndex) &&
                              nonlinear_it_index <
                                MatrixFreePDE<dim, degree>::userInputs
                                  .nonlinear_solver_parameters.getMaxIterations())
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
                          if (this->fields[fieldIndex].type == SCALAR)
                            {
                              *this->solutionSet[fieldIndex] += this->dU_scalar;
                            }
                          else
                            {
                              *this->solutionSet[fieldIndex] += this->dU_vector;
                            }

                          if (this->currentIncrement % userInputs.skip_print_steps == 0)
                            {
                              double dU_norm;
                              if (this->fields[fieldIndex].type == SCALAR)
                                {
                                  dU_norm = this->dU_scalar.l2_norm();
                                }
                              else
                                {
                                  dU_norm = this->dU_vector.l2_norm();
                                }
                              snprintf(buffer,
                                       sizeof(buffer),
                                       "field '%2s' [linear solve]: initial "
                                       "residual:%12.6e, current residual:%12.6e, "
                                       "nsteps:%u, tolerance criterion:%12.6e, "
                                       "solution: %12.6e, dU: %12.6e\n",
                                       this->fields[fieldIndex].name.c_str(),
                                       this->residualSet[fieldIndex]->l2_norm(),
                                       solver_control.last_value(),
                                       solver_control.last_step(),
                                       solver_control.tolerance(),
                                       this->solutionSet[fieldIndex]->l2_norm(),
                                       dU_norm);
                              this->pcout << buffer;
                            }
                        }
                    }
                }
              else if (this->fields[fieldIndex].pdetype == AUXILIARY)
                {
                  if (userInputs.var_nonlinear[fieldIndex] || nonlinear_it_index == 0)
                    {
                      // If the equation for this field is nonlinear, save the
                      // old solution
                      if (userInputs.var_nonlinear[fieldIndex])
                        {
                          if (this->fields[fieldIndex].type == SCALAR)
                            {
                              this->dU_scalar = *this->solutionSet[fieldIndex];
                            }
                          else
                            {
                              this->dU_vector = *this->solutionSet[fieldIndex];
                            }
                        }

                      this->updateExplicitSolution(fieldIndex);

                      // Apply Boundary conditions
                      this->applyBCs(fieldIndex);

                      // Print update to screen
                      if (this->currentIncrement % userInputs.skip_print_steps == 0)
                        {
                          snprintf(buffer,
                                   sizeof(buffer),
                                   "field '%2s' [auxiliary solve]: current solution: "
                                   "%12.6e, current residual:%12.6e\n",
                                   this->fields[fieldIndex].name.c_str(),
                                   this->solutionSet[fieldIndex]->l2_norm(),
                                   this->residualSet[fieldIndex]->l2_norm());
                          this->pcout << buffer;
                        }

                      // Check to see if this individual variable has converged
                      if (userInputs.var_nonlinear[fieldIndex])
                        {
                          if (MatrixFreePDE<dim, degree>::userInputs
                                .nonlinear_solver_parameters.getToleranceType(
                                  fieldIndex) == ABSOLUTE_SOLUTION_CHANGE)
                            {
                              double diff;

                              if (this->fields[fieldIndex].type == SCALAR)
                                {
                                  this->dU_scalar -= *this->solutionSet[fieldIndex];
                                  diff = this->dU_scalar.l2_norm();
                                }
                              else
                                {
                                  this->dU_vector -= *this->solutionSet[fieldIndex];
                                  diff = this->dU_vector.l2_norm();
                                }
                              if (this->currentIncrement % userInputs.skip_print_steps ==
                                  0)
                                {
                                  this->pcout << "Relative difference between nonlinear "
                                                 "iterations: "
                                              << diff << " " << nonlinear_it_index << " "
                                              << this->currentIncrement << std::endl;
                                }

                              if (diff > MatrixFreePDE<dim, degree>::userInputs
                                           .nonlinear_solver_parameters.getToleranceValue(
                                             fieldIndex) &&
                                  nonlinear_it_index <
                                    MatrixFreePDE<dim, degree>::userInputs
                                      .nonlinear_solver_parameters.getMaxIterations())
                                {
                                  nonlinear_it_converged = false;
                                }
                            }
                          else
                            {
                              std::cerr << "PRISMS-PF Error: Nonlinear solver "
                                           "tolerance types other than ABSOLUTE_CHANGE "
                                           "have yet to be implemented."
                                        << std::endl;
                            }
                        }
                    }
                }

              // check if solution is nan
              if (!numbers::is_finite(this->solutionSet[fieldIndex]->l2_norm()))
                {
                  snprintf(buffer,
                           sizeof(buffer),
                           "ERROR: field '%s' solution is NAN. exiting.\n\n",
                           this->fields[fieldIndex].name.c_str());
                  this->pcout << buffer;
                  exit(-1);
                }
            }

          nonlinear_it_index++;
        }
    }

  if (this->currentIncrement % userInputs.skip_print_steps == 0)
    {
      this->pcout << "wall time: " << time.wall_time() << "s\n";
    }
  // log time
  this->computing_timer.leave_subsection("matrixFreePDE: solveIncrements");
}
