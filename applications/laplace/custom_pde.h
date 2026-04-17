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
  using ScalarGrad  = dealii::Tensor<1, dim, ScalarValue>;
  using ScalarHess  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;
  using VectorGrad  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorHess  = dealii::Tensor<3, dim, ScalarValue>;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  /**
   * @brief Constructor.
   */
  explicit CustomPDE(const UserInputParameters<dim> &_user_inputs,
                     PhaseFieldTools<dim>           &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
  {}

  void
  set_dirichlet([[maybe_unused]] const unsigned int       &index,
                [[maybe_unused]] const unsigned int       &boundary_id,
                [[maybe_unused]] const unsigned int       &component,
                [[maybe_unused]] const dealii::Point<dim> &point,
                [[maybe_unused]] number                   &scalar_value,
                [[maybe_unused]] number &vector_component_value) const override
  {
    [[maybe_unused]] const double x = (dim > 0) ? point[0] : 0.0;
    [[maybe_unused]] const double y = (dim > 1) ? point[1] : 0.0;
    [[maybe_unused]] const double z = (dim > 2) ? point[2] : 0.0;

    if (boundary_id == 0) // left
      {
        using std::cos;
        scalar_value = 50.0 * cos(x / 4.0) * cos(y / 4.0);
        return;
      }
    if (boundary_id == 1) // right
      {
        scalar_value = y / 2.0;
        return;
      }
    if (boundary_id == 2) // bottom
      {
        scalar_value = x * (100.0 - x) / 50.0;
        return;
      }
    if (boundary_id == 3) // top
      {
        scalar_value = 0.0;
        return;
      }
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 1) // linear rhs
      {
        variable_list.set_value_term(0, ScalarValue(0.0));
      }
  }

  void
  compute_lhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 1) // linear lhs
      {
        ScalarGrad ux = variable_list.template get_gradient<Scalar, LHS>(0);
        variable_list.set_gradient_term(0, ux);
      }
  }
};

PRISMS_PF_END_NAMESPACE
