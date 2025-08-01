// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/pde_operator.h>
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

  /**
   * @brief Constructor.
   */
  explicit CustomPDE(const UserInputParameters<dim> &_user_inputs)
    : PDEOperator<dim, degree, number>(_user_inputs)
  {}

private:
  /**
   * @brief User-implemented class for the initial conditions.
   */
  void
  set_initial_condition(const unsigned int       &index,
                        const unsigned int       &component,
                        const dealii::Point<dim> &point,
                        double                   &scalar_value,
                        double                   &vector_component_value) const override;

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
    Types::Index solve_block) const override;

  /**
   * @brief User-implemented class for the RHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_rhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    Types::Index                                               solve_block,
    Types::Index index = Numbers::invalid_index) const override;

  /**
   * @brief User-implemented class for the LHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_lhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    Types::Index                                               solve_block,
    Types::Index index = Numbers::invalid_index) const override;

  /**
   * @brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  void
  compute_postprocess_explicit_rhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    Types::Index solve_block) const override;

  number McV =
    this->get_user_inputs().get_user_constants().get_model_constant_double("McV");
  number Mn1V =
    this->get_user_inputs().get_user_constants().get_model_constant_double("Mn1V");
  number Mn2V =
    this->get_user_inputs().get_user_constants().get_model_constant_double("Mn2V");
  number Mn3V =
    this->get_user_inputs().get_user_constants().get_model_constant_double("Mn3V");
  dealii::Tensor<2, dim, number> Kn1 =
    this->get_user_inputs().get_user_constants().get_model_constant_rank_2_tensor("Kn1");
  dealii::Tensor<2, dim, number> Kn2 =
    this->get_user_inputs().get_user_constants().get_model_constant_rank_2_tensor("Kn2");
  dealii::Tensor<2, dim, number> Kn3 =
    this->get_user_inputs().get_user_constants().get_model_constant_rank_2_tensor("Kn3");

  bool n_dependent_stiffness =
    this->get_user_inputs().get_user_constants().get_model_constant_bool(
      "n_dependent_stiffness");

  dealii::Tensor<2, dim, number> sfts_const1 =
    this->get_user_inputs().get_user_constants().get_model_constant_rank_2_tensor(
      "sfts_const1");
  dealii::Tensor<2, dim, number> sfts_const2 =
    this->get_user_inputs().get_user_constants().get_model_constant_rank_2_tensor(
      "sfts_const2");
  dealii::Tensor<2, dim, number> sfts_const3 =
    this->get_user_inputs().get_user_constants().get_model_constant_rank_2_tensor(
      "sfts_const3");

  number A4 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("A4");
  number A3 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("A3");
  number A2 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("A2");
  number A1 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("A1");
  number A0 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("A0");

  number B2 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("B2");
  number B1 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("B1");
  number B0 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("B0");

  dealii::Tensor<2, voigt_tensor_size<dim>, number> CIJ_Mg =
    this->get_user_inputs().get_user_constants().get_model_constant_elasticity_tensor(
      "CIJ_Mg");
  dealii::Tensor<2, voigt_tensor_size<dim>, number> CIJ_Beta =
    this->get_user_inputs().get_user_constants().get_model_constant_elasticity_tensor(
      "CIJ_Beta");
};

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
