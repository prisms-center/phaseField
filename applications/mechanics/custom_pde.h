// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>

#include <random>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public PDEOperatorBase<dim, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using ScalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using VectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  /**
   * @brief Constructor.
   */
  explicit CustomPDE(const UserInputParameters<dim> &_user_inputs)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , m_well(get_user_inputs().user_constants.get_model_constant_double("m_well"))
    , kappa(get_user_inputs().user_constants.get_model_constant_double("kappa"))
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const Point<dim>   &p,
                        [[maybe_unused]] const unsigned int &index,
                        [[maybe_unused]] const double       &scalar_IC,
                        [[maybe_unused]] Vector<double>     &vector_IC) override;
  {}

  void
  compute_rhs(FieldContainer<dim, degree, number>      &variable_list,
              const Point<dim, VectorizedArray<double>> q_point_loc,
              const VectorizedArray<double>             element_volume) const override;
  {
    if (solve_group_id == 1) // explicit
      {
        VectorGrad              ux = variable_list.template get_vector_gradient(0);
        VectorGrad              eqx_u;
        VectorizedArray<double> E[dim][dim], S[dim][dim];
        for (unsigned int i = 0; i < dim; i++)
          {
            for (unsigned int j = 0; j < dim; j++)
              {
                E[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]);
              }
          }
        computeStress<dim>(CIJ, E, S);

        for (unsigned int i = 0; i < dim; i++)
          {
            for (unsigned int j = 0; j < dim; j++)
              {
                eqx_u[i][j] = -S[i][j];
              }
          }
        variable_list.set_vector_gradient_term_RHS(0, eqx_u);
      }
  }

  number m_well;
  number kappa;
};

PRISMS_PF_END_NAMESPACE
