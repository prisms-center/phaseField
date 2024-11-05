#include "../../include/SolverParameters.h"

SolverToleranceType
SolverParametersBase::getToleranceType(unsigned int index)
{
  return tolerance_type_list.at(getEquationIndex(index));
}

double
SolverParametersBase::getToleranceValue(unsigned int index)
{
  return tolerance_value_list.at(getEquationIndex(index));
}

unsigned int
SolverParametersBase::getEquationIndex(unsigned int global_index)
{
  for (unsigned int i = 0; i < var_index_list.size(); i++)
    {
      if (var_index_list.at(i) == global_index)
        {
          return i;
        }
    }
  std::cerr << "PRISMS-PF Error: Attempted access of a parameter for the "
               "nonlinear solver for an ineligible variable index.\n";
  abort();
}

void
LinearSolverParameters::loadParameters(unsigned int        _var_index,
                                       SolverToleranceType _tolerance_type,
                                       double              _tolerance_value,
                                       unsigned int        _max_iterations)
{
  var_index_list.push_back(_var_index);
  tolerance_type_list.push_back(_tolerance_type);
  tolerance_value_list.push_back(_tolerance_value);
  max_iterations_list.push_back(_max_iterations);
}

unsigned int
LinearSolverParameters::getMaxIterations(unsigned int index)
{
  return max_iterations_list.at(getEquationIndex(index));
}

void
NonlinearSolverParameters::loadParameters(unsigned int        _var_index,
                                          SolverToleranceType _tolerance_type,
                                          double              _tolerance_value,
                                          bool                _backtrack_damping_flag,
                                          double              _backtrack_step_modifier,
                                          double _backtrack_residual_decrease_coeff,
                                          double _default_dampling_coefficient,
                                          bool   _laplace_for_initial_guess)
{
  var_index_list.push_back(_var_index);
  tolerance_type_list.push_back(_tolerance_type);
  tolerance_value_list.push_back(_tolerance_value);
  backtrack_damping_flag_list.push_back(_backtrack_damping_flag);
  backtrack_step_modifier_list.push_back(_backtrack_step_modifier);
  backtrack_residual_decrease_coeff_list.push_back(_backtrack_residual_decrease_coeff);
  default_damping_coefficient_list.push_back(_default_dampling_coefficient);
  laplace_for_initial_guess_list.push_back(_laplace_for_initial_guess);
}

bool
NonlinearSolverParameters::getBacktrackDampingFlag(unsigned int index)
{
  return backtrack_damping_flag_list.at(getEquationIndex(index));
}

double
NonlinearSolverParameters::getBacktrackStepModifier(unsigned int index)
{
  return backtrack_step_modifier_list.at(getEquationIndex(index));
}

double
NonlinearSolverParameters::getBacktrackResidualDecreaseCoeff(unsigned int index)
{
  return backtrack_residual_decrease_coeff_list.at(getEquationIndex(index));
}

double
NonlinearSolverParameters::getDefaultDampingCoefficient(unsigned int index)
{
  return default_damping_coefficient_list.at(getEquationIndex(index));
}

bool
NonlinearSolverParameters::getLaplaceInitializationFlag(unsigned int index)
{
  return laplace_for_initial_guess_list.at(getEquationIndex(index));
}

void
NonlinearSolverParameters::setMaxIterations(unsigned int _max_iterations)
{
  max_iterations = _max_iterations;
}

double
NonlinearSolverParameters::getMaxIterations() const
{
  return max_iterations;
}
