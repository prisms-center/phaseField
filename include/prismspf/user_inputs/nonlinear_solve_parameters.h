#ifndef nonlinear_solve_parameters_h
#define nonlinear_solve_parameters_h

#include <prismspf/config.h>
#include <prismspf/user_inputs/linear_solve_parameters.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Class that stores relevant nonlinear solve information of a certain field
 */
class nonlinearSolverParameters
{
public:
  // Nonlinear step length
  mutable double step_length = 1.0;

  // Max number of iterations for the nonlinear solve
  unsigned int max_iterations = 100;
};

/**
 * \brief Class that holds nonlinear solver parameters.
 */
struct nonlinearSolveParameters
{
public:
  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Map of nonlinear solve parameters for fields that require them
  std::map<unsigned int, nonlinearSolverParameters> nonlinear_solve;
};

inline void
nonlinearSolveParameters::print_parameter_summary() const
{
  if (!nonlinear_solve.empty())
    {
      prisms::conditionalOStreams::pout_summary()
        << "================================================\n"
        << "  Nonlinear Solve Parameters\n"
        << "================================================\n";

      for (const auto &[index, nonlinear_solver_parameters] : nonlinear_solve)
        {
          prisms::conditionalOStreams::pout_summary()
            << "Index: " << index << "\n"
            << "  Max iterations: " << nonlinear_solver_parameters.max_iterations << "\n"
            << "  Step length: " << nonlinear_solver_parameters.step_length << "\n";
        }

      prisms::conditionalOStreams::pout_summary() << "\n" << std::flush;
    }
}

PRISMS_PF_END_NAMESPACE

#endif