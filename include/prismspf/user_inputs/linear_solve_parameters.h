#ifndef linear_solve_parameters_h
#define linear_solve_parameters_h

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>

#include <map>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Class that stores relevant linear solve information of a certain field
 */
struct linearSolverParameters
{
public:
  // Solver tolerance
  double tolerance = 1.0e-6;

  // Solver tolerance type
  solverToleranceType tolerance_type = solverToleranceType::RELATIVE_RESIDUAL_CHANGE;

  // Max number of iterations for the linear solve
  unsigned int max_iterations = 100;

  // Preconditioner
  preconditionerType preconditioner = preconditionerType::GMG;

  // Smoothing range for eigenvalues. This denotes the lower bound of eigenvalues that are
  // smoothed [1.2 λ^max / smoothing_range, 1.2 λ^max], where λ^max is the estimated
  // maximum eigenvalue. A choice between 5 and 20 is usually useful when the
  // preconditioner is used as a smoother in multigrid.
  double smoothing_range = 15.0;

  // Polynomial degree for the Chebyshev smoother
  unsigned int smoother_degree = 5;

  // Maximum number of CG iterations used to find the maximum eigenvalue
  unsigned int eig_cg_n_iterations = 10;
};

/**
 * \brief Class that holds linear solver parameters.
 */
struct linearSolveParameters
{
public:
  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Map of linear solve parameters for fields that require them
  std::map<unsigned int, linearSolverParameters> linear_solve;
};

inline void
linearSolveParameters::print_parameter_summary() const
{
  if (!linear_solve.empty())
    {
      prisms::conditionalOStreams::pout_summary()
        << "================================================\n"
        << "  Linear Solve Parameters\n"
        << "================================================\n";

      for (const auto &[index, linear_solver_parameters] : linear_solve)
        {
          prisms::conditionalOStreams::pout_summary()
            << "Index: " << index << "\n"
            << "  Tolerance: " << linear_solver_parameters.tolerance << "\n"
            << "  Type: " << to_string(linear_solver_parameters.tolerance_type) << "\n"
            << "  Max iterations: " << linear_solver_parameters.max_iterations << "\n"
            << "  Preconditioner: " << to_string(linear_solver_parameters.preconditioner)
            << "\n";

          if (linear_solver_parameters.preconditioner == preconditionerType::GMG)
            {
              prisms::conditionalOStreams::pout_summary()
                << "  Smoothing range: " << linear_solver_parameters.smoothing_range
                << "\n"
                << "  Smoother degree: " << linear_solver_parameters.smoother_degree
                << "\n"
                << "  Max eigenvalue CG iterations: "
                << linear_solver_parameters.eig_cg_n_iterations << "\n";
            }
        }

      prisms::conditionalOStreams::pout_summary() << "\n" << std::flush;
    }
}

PRISMS_PF_END_NAMESPACE

#endif