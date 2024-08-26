// Definition of the variables in the model
#define num_var 1
#define variable_name \
  {                   \
    "n"               \
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

// define Allen-Cahn parameters
#define MnV 1.0
#define KnV 4.0

// define the free energy and its derivative with respect to n
#define fV (n * n * n * n - 2.0 * n * n * n + n * n)
#define fnV (4.0 * n * (n - 1.0) * (n - 0.5))

// Define required residuals
#define rnV (n - constV(timeStep * MnV) * fnV)
#define rnxV (constV(-timeStep * KnV * MnV) * nx)

// ---------------------------------------------

template <int dim>
void
generalizedProblem<dim>::residualRHS(
  const std::vector<modelVariable<dim>>              &modelVariablesList,
  std::vector<modelResidual<dim>>                    &modelResidualsList,
  dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const
{
  // n
  scalarvalueType n  = modelVariablesList[0].scalarValue;
  scalargradType  nx = modelVariablesList[0].scalarGrad;

  modelResidualsList[0].scalarValueResidual = rnV;
  modelResidualsList[0].scalarGradResidual  = rnxV;
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
{
  scalarvalueType total_energy_density = constV(0.0);

  // n
  scalarvalueType n  = modelVarList[0].scalarValue;
  scalargradType  nx = modelVarList[0].scalarGrad;

  scalarvalueType f_chem = fV;

  scalarvalueType f_grad = constV(0.0);

  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += constV(0.5 * KnV) * nx[i] * nx[j];
        }
    }

  total_energy_density = f_chem + f_grad;

  assembler_lock.acquire();
  for (unsigned i = 0; i < n.size(); i++)
    {
      // For some reason, some of the values in this loop
      if (n[i] > 1.0e-10)
        {
          this->energy += total_energy_density[i] * JxW_value[i];
          this->energy_components[0] += f_chem[i] * JxW_value[i];
          this->energy_components[1] += f_grad[i] * JxW_value[i];
        }
    }
  assembler_lock.release();
}
