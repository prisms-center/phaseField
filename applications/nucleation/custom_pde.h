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
    Types::Index current_index = Numbers::invalid_index) const override;

  /**
   * @brief User-implemented class for the LHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_lhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    const dealii::VectorizedArray<number>                     &element_volume,
    Types::Index                                               solve_block,
    Types::Index current_index = Numbers::invalid_index) const override;

  /**
   * @brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  void
  compute_postprocess_explicit_rhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    const dealii::VectorizedArray<number>                     &element_volume,
    Types::Index solve_block) const override;

  void
  seed_nucleus(const dealii::Point<dim, ScalarValue> &q_point_loc,
               ScalarValue                           &source_term,
               ScalarValue                           &gamma) const;

  double c_avg =
    this->get_user_inputs().get_user_constants().get_model_constant_double("c_avg");
  double McV =
    this->get_user_inputs().get_user_constants().get_model_constant_double("McV");
  double MnV =
    this->get_user_inputs().get_user_constants().get_model_constant_double("MnV");
  double KnV =
    this->get_user_inputs().get_user_constants().get_model_constant_double("KnV");
  double W_barrier =
    this->get_user_inputs().get_user_constants().get_model_constant_double("W_barrier");
  double A0 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("A0");
  double A2 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("A2");
  double calmin =
    this->get_user_inputs().get_user_constants().get_model_constant_double("calmin");
  double B0 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("B0");
  double B2 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("B2");
  double cbtmin =
    this->get_user_inputs().get_user_constants().get_model_constant_double("cbtmin");

  double k1 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("k1");
  double k2 =
    this->get_user_inputs().get_user_constants().get_model_constant_double("k2");
  double tau =
    this->get_user_inputs().get_user_constants().get_model_constant_double("tau");
  double epsilon =
    this->get_user_inputs().get_user_constants().get_model_constant_double("epsilon");
  double r_nuc =
    this->get_user_inputs().get_user_constants().get_model_constant_double("r_nuc");
  double r_freeze =
    this->get_user_inputs().get_user_constants().get_model_constant_double("r_freeze");
  double seeding_duration =
    this->get_user_inputs().get_user_constants().get_model_constant_double(
      "seeding_duration");

  // Interface coefficient
  double interface_coeff = std::sqrt(2.0 * KnV / W_barrier);
};

PRISMS_PF_END_NAMESPACE