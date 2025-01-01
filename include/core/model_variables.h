#ifndef INCLUDE_MODELVARIABLE_H_
#define INCLUDE_MODELVARIABLE_H_

#include <deal.II/matrix_free/evaluation_flags.h>

struct variable_info
{
  using EvalFlags = dealii::EvaluationFlags::EvaluationFlags;

  unsigned int global_var_index = 0;
  bool         is_scalar        = true;
  bool         var_needed       = false;

  EvalFlags evaluation_flags = dealii::EvaluationFlags::nothing;
  EvalFlags residual_flags   = dealii::EvaluationFlags::nothing;
};

#endif
