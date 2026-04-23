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
    , stiffness(get_user_inputs().user_constants.get_elasticity_tensor("stiffness"))
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {}

  void
  set_dirichlet([[maybe_unused]] const unsigned int       &index,
                [[maybe_unused]] const unsigned int       &boundary_id,
                [[maybe_unused]] const unsigned int       &component,
                [[maybe_unused]] const dealii::Point<dim> &point,
                [[maybe_unused]] number                   &scalar_value,
                [[maybe_unused]] number &vector_component_value) const override
  {
    scalar_value           = 0.0;
    vector_component_value = 0.0;
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 1) // linear rhs
      {
        ScalarValue dist_from_inclusion = variable_list.get_q_point_location().norm();
        ScalarValue inclusion_radius(10.0);

        VectorGrad  transformation_strain;
        ScalarValue strain_value =
          0.01 * (0.5 + 0.5 * std::tanh(10.0 * (dist_from_inclusion - inclusion_radius)));
        for (unsigned int i = 0; i < dim; i++)
          {
            transformation_strain[i][i] = strain_value;
          }
        VectorGrad stress;
        compute_stress<dim, ScalarValue>(stiffness, transformation_strain, stress);
        variable_list.set_gradient_term(0, -stress);
      }
  }

  void
  compute_lhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 1) // linear lhs
      {
        VectorGrad ux = variable_list.template get_symmetric_gradient<Vector, LHS>(0);
        VectorGrad stress;
        compute_stress<dim, ScalarValue>(stiffness, ux, stress);
        variable_list.set_gradient_term(0, stress);
      }
  }

  constexpr static unsigned int CIJ_tensor_size = (2 * dim) - 1 + (dim / 3);

  dealii::Tensor<2, CIJ_tensor_size, number> stiffness;
};

PRISMS_PF_END_NAMESPACE
