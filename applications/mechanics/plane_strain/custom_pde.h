// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/utilities/mechanics.h>

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
  CustomPDE(const UserInputParameters<dim> &_user_inputs, PhaseFieldTools<dim> &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , stiffness(
        get_user_inputs().user_constants.get_elasticity_tensor_plane_strain("stiffness"))
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
                [[maybe_unused]] const SimulationTimer    &sim_timer,
                [[maybe_unused]] number                   &scalar_value,
                [[maybe_unused]] number &vector_component_value) const override
  {
    // Zero except for x-component of u on x=0 face
    if (index == 0 && boundary_id == RectangularMesh<dim>::Boundary::Left &&
        component == 0)
      {
        vector_component_value = -1.0;
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
        variable_list.set_value_term(0, VectorValue());
      }
    else if (solve_block_id == 2) // post-processing
      {
        VectorGrad strain =
          variable_list.template get_symmetric_gradient<Vector, Current>(0);
        VectorGrad  stress {};
        ScalarValue stress_zz {};
        Mechanics::compute_stress<dim, StressState::PlaneStrain, ScalarValue>(stiffness,
                                                                              strain,
                                                                              stress,
                                                                              stress_zz);

        variable_list.set_value_term(1, strain[0][0]);
        variable_list.set_value_term(2, strain[1][1]);
        variable_list.set_value_term(3, 2.0 * strain[0][1]);
        variable_list.set_value_term(4, stress[0][0]);
        variable_list.set_value_term(5, stress[1][1]);
        variable_list.set_value_term(6, stress[0][1]);
        variable_list.set_value_term(7, stress_zz);
      }
  }

  void
  compute_lhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 1) // linear lhs
      {
        VectorGrad strain = variable_list.template get_symmetric_gradient<Vector, LHS>(0);
        VectorGrad stress {};
        ScalarValue stress_zz {};
        // if strain_zz = 0
        Mechanics::compute_stress<dim, StressState::PlaneStrain, ScalarValue>(stiffness,
                                                                              strain,
                                                                              stress,
                                                                              stress_zz);
        // if strain_zz is not 0
        //  Mechanics::compute_stress<dim, StressState::PlaneStrain,
        //  ScalarValue>(stiffness, strain, strain_zz, stress, stress_zz);

        variable_list.set_gradient_term(0, stress);
      }
  }

  dealii::Tensor<2, Mechanics::voigt_size<dim, StressState::PlaneStrain>, number>
    stiffness;
};

PRISMS_PF_END_NAMESPACE
