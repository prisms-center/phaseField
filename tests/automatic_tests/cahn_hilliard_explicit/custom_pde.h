// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/pde_operator.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This is a derived class of `MatrixFreeOperator` where the user implements their
 * PDEs.
 *
 * \tparam dim The number of dimensions in the problem.
 * \tparam degree The polynomial degree of the shape functions.
 * \tparam number Datatype to use. Either double or float.
 */
template <unsigned int dim, unsigned int degree, typename number>
class customPDE : public PDEOperator<dim, degree, number>
{
public:
  using scalarValue = dealii::VectorizedArray<number>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using scalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using vectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;

  /**
   * \brief Constructor.
   */
  explicit customPDE(const UserInputParameters<dim> &_user_inputs)
    : PDEOperator<dim, degree, number>(_user_inputs)
  {}

private:
  /**
   * \brief User-implemented class for the initial conditions.
   */
  void
  set_initial_condition(const unsigned int       &index,
                        const unsigned int       &component,
                        const dealii::Point<dim> &point,
                        double                   &scalar_value,
                        double                   &vector_component_value) const override;

  /**
   * \brief User-implemented class for nonuniform boundary conditions.
   */
  virtual void
  set_nonuniform_dirichlet(const unsigned int       &index,
                           const unsigned int       &boundary_id,
                           const unsigned int       &component,
                           const dealii::Point<dim> &point,
                           number                   &scalar_value,
                           number &vector_component_value) const override;

  /**
   * \brief User-implemented class for the RHS of explicit equations.
   */
  void
  compute_explicit_rhs(VariableContainer<dim, degree, number> &variable_list,
                       const dealii::Point<dim, dealii::VectorizedArray<number>>
                         &q_point_loc) const override;

  /**
   * \brief User-implemented class for the RHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_rhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    Types::Index current_index = Numbers::invalid_index) const override;

  /**
   * \brief User-implemented class for the LHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_lhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    Types::Index current_index = Numbers::invalid_index) const override;

  /**
   * \brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  void
  compute_postprocess_explicit_rhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
    const override;

  number McV =
    this->get_user_inputs().get_user_constants().get_model_constant_double("McV");
  number KcV =
    this->get_user_inputs().get_user_constants().get_model_constant_double("KcV");
};

PRISMS_PF_END_NAMESPACE