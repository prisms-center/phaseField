// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This is a derived class of `matrixFreeOperator` where the user implements their
 * PDEs.
 *
 * \tparam dim The number of dimensions in the problem.
 * \tparam degree The polynomial degree of the shape functions.
 * \tparam number Datatype to use. Either double or float.
 */
template <int dim, int degree, typename number>
class customPDE : public matrixFreeOperator<dim, degree, number>
{
public:
  using scalarValue = dealii::VectorizedArray<number>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using scalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using vectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;

  /**
   * \brief Constructor for concurrent solves.
   */
  customPDE(const userInputParameters<dim>                   &_user_inputs,
            const std::map<unsigned int, variableAttributes> &subset_attributes)
    : matrixFreeOperator<dim, degree, number>(_user_inputs, subset_attributes)
  {}

  /**
   * \brief Constructor for single solves.
   */
  customPDE(const userInputParameters<dim>                   &_user_inputs,
            const unsigned int                               &_current_index,
            const std::map<unsigned int, variableAttributes> &subset_attributes)
    : matrixFreeOperator<dim, degree, number>(_user_inputs,
                                              _current_index,
                                              subset_attributes)
  {}

private:
  /**
   * \brief User-implemented class for the RHS of explicit equations.
   */
  void
  compute_explicit_RHS(variableContainer<dim, degree, number> &variable_list,
                       const dealii::Point<dim, dealii::VectorizedArray<number>>
                         &q_point_loc) const override;

  /**
   * \brief User-implemented class for the RHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_RHS(variableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, dealii::VectorizedArray<number>>
                            &q_point_loc) const override;

  /**
   * \brief User-implemented class for the LHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_LHS(variableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, dealii::VectorizedArray<number>>
                            &q_point_loc) const override;

  /**
   * \brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  void
  compute_postprocess_explicit_RHS(
    variableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
    const override;

  number McV  = this->user_inputs.user_constants.get_model_constant_double("McV");
  number Mn1V = this->user_inputs.user_constants.get_model_constant_double("Mn1V");
  dealii::Tensor<2, dim, number> Kn1 =
    this->user_inputs.user_constants.get_model_constant_rank_2_tensor("Kn1");
  number W = this->user_inputs.user_constants.get_model_constant_double("W");
  bool   n_dependent_stiffness =
    this->user_inputs.user_constants.get_model_constant_bool("n_dependent_stiffness");
  dealii::Tensor<2, dim, number> sfts_linear1 =
    this->user_inputs.user_constants.get_model_constant_rank_2_tensor("sfts_linear1");
  dealii::Tensor<2, dim, number> sfts_const1 =
    this->user_inputs.user_constants.get_model_constant_rank_2_tensor("sfts_const1");

  number A2 = this->user_inputs.user_constants.get_model_constant_double("A2");
  number A1 = this->user_inputs.user_constants.get_model_constant_double("A1");
  number A0 = this->user_inputs.user_constants.get_model_constant_double("A0");
  number B2 = this->user_inputs.user_constants.get_model_constant_double("B2");
  number B1 = this->user_inputs.user_constants.get_model_constant_double("B1");
  number B0 = this->user_inputs.user_constants.get_model_constant_double("B0");

  dealii::Tensor<2, (2 * dim) - 1 + (dim / 3), number> CIJ_Mg =
    this->user_inputs.user_constants.get_model_constant_elasticity_tensor("CIJ_Mg");
  dealii::Tensor<2, (2 * dim) - 1 + (dim / 3), number> CIJ_Beta =
    this->user_inputs.user_constants.get_model_constant_elasticity_tensor("CIJ_Beta");

  bool c_dependent_misfit;
};

PRISMS_PF_END_NAMESPACE