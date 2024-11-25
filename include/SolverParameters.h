#ifndef INCLUDE_NONLINEARSOLVERPARAMETERS_H_
#define INCLUDE_NONLINEARSOLVERPARAMETERS_H_

#include <iostream>
#include <vector>

enum SolverToleranceType
{
  ABSOLUTE_RESIDUAL,
  RELATIVE_RESIDUAL_CHANGE,
  ABSOLUTE_SOLUTION_CHANGE
};

/**
 * This is a base class that holds parameters related to a numerical solver
 * (either linear or nonlinear) As in many other PRISMS-PF classes, access to
 * the parameters for the separate governing equations is specified via the
 * global variable index. Requests for solver parameters for inappropriate
 * governing equations return errors.
 */
class SolverParametersBase
{
public:
  /**
   * Method to get the tolerance type for one of the governing equations.
   */
  SolverToleranceType
  getToleranceType(unsigned int index);

  /**
   * Method to get the tolerance value for one of the governing equations.
   */
  double
  getToleranceValue(unsigned int index);

protected:
  std::vector<unsigned int> var_index_list;

  std::vector<SolverToleranceType> tolerance_type_list;
  std::vector<double>              tolerance_value_list;

  /**
   * Method to get the internal nonlinearSolverParameters index from the global
   * variable index.
   */
  unsigned int
  getEquationIndex(unsigned int global_index);
};

/**
 * This class holds all of the parameters for the linear solver. Requests for
 * linear solver parameters for governing equations without a matrix inversion
 * return errors.
 */
class LinearSolverParameters : public SolverParametersBase
{
public:
  /**
   * Method to load the parameters for one governing equation into the class.
   * This must be defined for all classes derived from SolverParametersBase.
   */
  void
  loadParameters(unsigned int        _var_index,
                 SolverToleranceType _tolerance_type,
                 double              _tolerance_value,
                 unsigned int        _max_iterations);

  /**
   * Method to get the maximum number of allowed iterations for the linear
   * solver.
   */
  unsigned int
  getMaxIterations(unsigned int index);

protected:
  std::vector<unsigned int> max_iterations_list;
};

/**
 * This class holds all of the parameters for the nonlinear solver. Requests for
 * nonlinear solver parameters for non-nonlinear governing equations return
 * errors.
 */
class NonlinearSolverParameters : public SolverParametersBase
{
public:
  /**
   * Method to get the maximum number of allowed iterations for the nonlinear
   * solver.
   */
  void
  setMaxIterations(unsigned int _max_iterations);

  /**
   * Method to get the maximum number of allowed iterations for the nonlinear
   * solver.
   */
  [[nodiscard]] double
  getMaxIterations() const;

  /**
   * Method to load the parameters for one governing equation into the class.
   * This must be defined for all classes derived from SolverParametersBase.
   */
  void
  loadParameters(unsigned int        _var_index,
                 SolverToleranceType _tolerance_type,
                 double              _tolerance_value,
                 bool                _backtrack_damping_flag,
                 double              _backtrack_step_modifier,
                 double              _backtrack_residual_decrease_coeff,
                 double              _default_dampling_coefficient,
                 bool                _laplace_for_initial_guess);

  /**
   * Method to get the backtrack line-search damping flag for one of the
   * governing equations.
   */
  bool
  getBacktrackDampingFlag(unsigned int index);

  /**
   * Method to get the backtrack line-search damping step size modifier for one
   * of the governing equations.
   */
  double
  getBacktrackStepModifier(unsigned int index);

  /**
   * Method to get the backtrack line-search damping residual decrease
   * coefficient for one of the governing equations.
   */
  double
  getBacktrackResidualDecreaseCoeff(unsigned int index);

  /**
   * Method to get the default damping coefficient for one of the governing
   * equations.
   */
  double
  getDefaultDampingCoefficient(unsigned int index);

  /**
   * Method to get the flag determining if the solution to Laplace's equation
   * should be used as the initial guess for one of the governing equations.
   */
  bool
  getLaplaceInitializationFlag(unsigned int index);

private:
  unsigned int max_iterations;

  std::vector<bool>   backtrack_damping_flag_list;
  std::vector<double> backtrack_step_modifier_list;
  std::vector<double> backtrack_residual_decrease_coeff_list;
  std::vector<double> default_damping_coefficient_list;
  std::vector<bool>   laplace_for_initial_guess_list;
};

#endif
