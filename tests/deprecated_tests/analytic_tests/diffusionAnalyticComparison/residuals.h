// Definition of the variables in the model
#define num_var 1
#define variable_name \
  {                   \
    "c"               \
  }
#define variable_type \
  {                   \
    "SCALAR"          \
  }
#define variable_eq_type \
  {                      \
    "PARABOLIC"          \
  }
#define need_val \
  {              \
    true         \
  }
#define need_grad \
  {               \
    true          \
  }
#define need_hess \
  {               \
    false         \
  }
#define need_val_residual \
  {                       \
    true                  \
  }
#define need_grad_residual \
  {                        \
    true                   \
  }

// define Fickian diffusion parameters
#define DcV 1.0

// define required residuals
#define rcV (c)
#define rcxV (constV(-DcV * timeStep) * cx)

// ---------------------------------------------

template <int dim>
void
generalizedProblem<dim>::residualRHS(
  const std::vector<modelVariable<dim>>              &modelVariablesList,
  std::vector<modelResidual<dim>>                    &modelResidualsList,
  dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const
{
  // c
  scalarvalueType c  = modelVariablesList[0].scalarValue;
  scalargradType  cx = modelVariablesList[0].scalarGrad;

  modelResidualsList[0].scalarValueResidual = rcV;
  modelResidualsList[0].scalarGradResidual  = rcxV;
}

template <int dim>
void
generalizedProblem<dim>::residualLHS(
  const std::vector<modelVariable<dim>>              &modelVarList,
  modelResidual<dim>                                 &modelRes,
  dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const
{}

template <int dim>
void
generalizedProblem<dim>::energyDensity(
  const std::vector<modelVariable<dim>>              &modelVarList,
  const dealii::VectorizedArray<double>              &JxW_value,
  dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc)
{}
