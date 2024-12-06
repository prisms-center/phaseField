// Model Variables Class

#ifndef INCLUDE_MODELVARIABLE_H_
#define INCLUDE_MODELVARIABLE_H_

#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/evaluation_flags.h>

struct variable_info
{
  unsigned int global_var_index = 0;
  bool         is_scalar        = true;
  bool         var_needed       = false;

  dealii::EvaluationFlags::EvaluationFlags evaluation_flags =
    dealii::EvaluationFlags::nothing;
  dealii::EvaluationFlags::EvaluationFlags residual_flags =
    dealii::EvaluationFlags::nothing;
};

#endif /* INCLUDE_MODELVARIABLE_H_ */
