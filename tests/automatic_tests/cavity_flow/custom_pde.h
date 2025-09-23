// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/pde_operator.h>
#include <prismspf/core/phase_field_tools.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This is a derived class of `VariableAttributeLoader` where the user implements
 * their variable attributes and field declarations.
 */
class CustomAttributeLoader : public VariableAttributeLoader
{
public:
  /**
   * @brief Destructor.
   */
  ~CustomAttributeLoader() override = default;

  /**
   * @brief User-implemented method where the variable attributes are set for all fields.
   */
  void
  load_variable_attributes() override;
};

/**
 * @brief This is a derived class of `MatrixFreeOperator` where the user implements their
 * PDEs.
 *
 * @tparam dim The number of dimensions in the problem.
 * @tparam degree The polynomial degree of the shape functions.
 * @tparam number Datatype to use. Either double or float.
 */
template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public PDEOperator<dim, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using ScalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using VectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;
  using PDEOperator<dim, degree, number>::get_user_inputs;
  using PDEOperator<dim, degree, number>::get_pf_tools;
  using PDEOperator<dim, degree, number>::get_timestep;

  /**
   * @brief Constructor.
   */
  explicit CustomPDE(const UserInputParameters<dim> &_user_inputs,
                     PhaseFieldTools<dim>           &_pf_tools)
    : PDEOperator<dim, degree, number>(_user_inputs, _pf_tools)
  {}

private:
  /**
   * @brief User-implemented class for the initial conditions.
   */
  void
  set_initial_condition(const unsigned int       &index,
                        const unsigned int       &component,
                        const dealii::Point<dim> &point,
                        number                   &scalar_value,
                        number                   &vector_component_value) const override;

  /**
   * @brief User-implemented class for nonuniform boundary conditions.
   */
  void
  set_nonuniform_dirichlet(const unsigned int       &index,
                           const unsigned int       &boundary_id,
                           const unsigned int       &component,
                           const dealii::Point<dim> &point,
                           number                   &scalar_value,
                           number &vector_component_value) const override;

  /**
   * @brief User-implemented class for the RHS of explicit equations.
   */
  void
  compute_explicit_rhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    const dealii::VectorizedArray<number>                     &element_volume,
    Types::Index solve_block) const override;

  /**
   * @brief User-implemented class for the RHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_rhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    const dealii::VectorizedArray<number>                     &element_volume,
    Types::Index                                               solve_block,
    Types::Index index = Numbers::invalid_index) const override;

  /**
   * @brief User-implemented class for the LHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_lhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    const dealii::VectorizedArray<number>                     &element_volume,
    Types::Index                                               solve_block,
    Types::Index index = Numbers::invalid_index) const override;

  /**
   * @brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  void
  compute_postprocess_explicit_rhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    const dealii::VectorizedArray<number>                     &element_volume,
    Types::Index solve_block) const override;

  const ScalarValue Re =
    get_user_inputs().get_user_constants().get_model_constant_double("Re");

  // Calculate the dynamic viscosity based on the user-inputted Reynolds number
  const ScalarValue nu = 1.0 / Re;
  const ScalarValue dt = get_timestep();

  // Values for stabilization term
  const ScalarValue size_modifier = std::sqrt(4.0 / M_PI) / degree;
  const ScalarValue time_contribution =
    dealii::Utilities::fixed_power<2>(1.0 / get_timestep());

  /**
   * @brief Compute the stabilization paramters
   */
  ScalarValue
  compute_stabilization_parameter(const VectorValue &velocity,
                                  const ScalarValue &element_volume) const
  {
    // Norm of the local velocity
    ScalarValue u_l2norm = 1.0e-12 + velocity.norm_square();

    // Stabilization parameter
    ScalarValue h = std::sqrt(element_volume) * size_modifier;

    ScalarValue stabilization_parameter =
      1.0 / std::sqrt(time_contribution + 4.0 * u_l2norm / h / h +
                      9.0 * dealii::Utilities::fixed_power<2>(4.0 * nu / h / h));

    return stabilization_parameter;
  }
};

PRISMS_PF_END_NAMESPACE
