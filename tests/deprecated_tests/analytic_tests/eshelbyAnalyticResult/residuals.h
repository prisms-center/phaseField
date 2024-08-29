// Define required residuals
// Definition of the variables in the model
#define num_var 1
#define variable_name \
  {                   \
    "u"               \
  }
#define variable_type \
  {                   \
    "VECTOR"          \
  }
#define variable_eq_type \
  {                      \
    "PARABOLIC"          \
  }
#define need_val \
  {              \
    false        \
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
    false                 \
  }
#define need_grad_residual \
  {                        \
    true                   \
  }

// Define Mechanical properties
// Mechanical symmetry of the material and stiffness parameters
#define MaterialModels \
  {                    \
    "ISOTROPIC"        \
  }
#define MaterialConstants \
  {                       \
    {                     \
      22.5, 0.3           \
    }                     \
  }

// ---------------------------------------------

template <int dim>
void
generalizedProblem<dim>::residualRHS(
  const std::vector<modelVariable<dim>>              &modelVariablesList,
  std::vector<modelResidual<dim>>                    &modelResidualsList,
  dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const
{
  // u
  vectorgradType ux = modelVariablesList[0].vectorGrad;
  vectorgradType Rux;

  dealii::VectorizedArray<double> sfts[dim][dim];

  dealii::VectorizedArray<double> dist;

  dist = std::sqrt(
    (q_point_loc[0] - constV(spanX / 2.0)) * (q_point_loc[0] - constV(spanX / 2.0)) +
    (q_point_loc[1] - constV(spanY / 2.0)) * (q_point_loc[1] - constV(spanY / 2.0)) +
    (q_point_loc[2] - constV(spanY / 2.0)) * (q_point_loc[2] - constV(spanY / 2.0)));

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          if (i == j)
            {
              sfts[i][j] =
                0.01 *
                (0.5 + 0.5 * (constV(1.0) - std::exp(-20.0 * (dist - constV(5.0)))) /
                         (constV(1.0) + std::exp(-20.0 * (dist - constV(5.0)))));
            }
          else
            {
              sfts[i][j] = 0.0;
            }
        }
    }

  // compute strain tensor
  dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]) - sfts[i][j];
        }
    }

  // compute stress tensor
  computeStress<dim>(CIJ_list[0], E, S);

  // compute residual
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          Rux[i][j] = -S[i][j];
        }
    }

  modelResidualsList[0].vectorGradResidual = Rux;
}

template <int dim>
void
generalizedProblem<dim>::residualLHS(
  const std::vector<modelVariable<dim>>              &modelVarList,
  modelResidual<dim>                                 &modelRes,
  dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const
{
  // u
  vectorgradType ux = modelVarList[0].vectorGrad;
  vectorgradType Rux;

  // compute strain tensor
  dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]);
        }
    }

  // compute stress tensor
  computeStress<dim>(CIJ_list[0], E, S);

  // compute residual
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          Rux[i][j] = S[i][j];
        }
    }

  modelRes.vectorGradResidual = Rux;
}

template <int dim>
void
generalizedProblem<dim>::energyDensity(
  const std::vector<modelVariable<dim>>              &modelVarList,
  const dealii::VectorizedArray<double>              &JxW_value,
  dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc)
{
  // u
  vectorgradType ux = modelVarList[0].vectorGrad;

  dealii::VectorizedArray<double> sfts[dim][dim];

  dealii::VectorizedArray<double> dist;

  dist = std::sqrt(
    (q_point_loc[0] - constV(spanX / 2.0)) * (q_point_loc[0] - constV(spanX / 2.0)) +
    (q_point_loc[1] - constV(spanY / 2.0)) * (q_point_loc[1] - constV(spanY / 2.0)) +
    (q_point_loc[2] - constV(spanY / 2.0)) * (q_point_loc[2] - constV(spanY / 2.0)));

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          if (i == j)
            {
              sfts[i][j] =
                0.01 *
                (0.5 + 0.5 * (constV(1.0) - std::exp(-20.0 * (dist - constV(10.0)))) /
                         (constV(1.0) + std::exp(-20.0 * (dist - constV(10.0)))));
            }
          else
            {
              sfts[i][j] = 0.0;
            }
        }
    }

  // compute strain tensor
  dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]) - sfts[i][j];
        }
    }

  // compute stress tensor
  computeStress<dim>(CIJ_list[0], E, S);

  scalarvalueType f_el = constV(0.0);

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_el += constV(0.5) * S[i][j] * E[i][j];
        }
    }

  assembler_lock.acquire();
  for (unsigned i = 0; i < f_el.size(); i++)
    {
      // For some reason, some of the values in this loop
      if (f_el[i] > 1.0e-10)
        {
          this->energy += f_el[i] * JxW_value[i];
        }
    }
  assembler_lock.release();
}
