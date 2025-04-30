// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

//#include <core/matrixFreePDE.h> //JM moved; maybe delete

//JM copied includes from precipitate app
#pragma once

#include <prismspf/core/pde_operator.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE //new namespace

//JM included for posterity
/**
 * \brief This is a derived class of `matrixFreeOperator` where the user implements their
 * PDEs.
 *
 * \tparam dim The number of dimensions in the problem.
 * \tparam degree The polynomial degree of the shape functions.
 * \tparam number Datatype to use. Either double or float.
 */

template <unsigned int dim, unsigned int degree, typename number> //JM unsigned ints and new typename
class customPDE : public PDEOperator<dim, degree, number>
{
public:
  //JM no namespace so declaring what to use explicitly
  using scalarValue = dealii::VectorizedArray<number>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using scalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using vectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;

  /**
  * \brief Constructor.
  */
  explicit customPDE(const userInputParameters<dim> &_user_inputs)
    : PDEOperator<dim, degree, number>(_user_inputs)
  {}

private:
//#include <core/typeDefs.h> //JM does not exist anymore
  /**
   * \brief User-implemented class for the RHS of explicit equations.
   */

  // Function to set the RHS of the governing equations for explicit time
  //JM new explicit RHS
  void
  compute_explicit_RHS(variableContainer<dim, degree, number> &variable_list,
                       const dealii::Point<dim, dealii::VectorizedArray<number>>
                         &q_point_loc) const override;

  /**
    * \brief User-implemented class for the RHS of nonexplicit equations.
    */
  //JM new nonexplicit RHS
  void
  compute_nonexplicit_RHS(
    variableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    types::index current_index = numbers::invalid_index) const override;
  
  /**
   * \brief User-implemented class for the LHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_LHS(
    variableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>>  &q_point_loc,
    types::index current_index = numbers::invalid_index) const override;


  /**
   * \brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  void
  compute_postprocess_explicit_RHS(
    variableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
    const override;

  constexpr static unsigned int  CIJ_tensor_size = 2 * dim - 1 + dim / 3; //JM modified to constexpr
  
  //JM syntax update for stiffness tensor and other constants
  dealii::Tensor<2, CIJ_tensor_size, number> CIJ = 
    this->get_user_inputs().user_constants.get_model_constant_elasticity_tensor("CIJ");

  number incRadius = this->get_user_inputs().user_constants.get_model_constant_double("incRadius");
  number poisson = this->get_user_inputs().user_constants.get_model_constant_double("poisson");
  dealii::Tensor<1, dim, number> center =
  this->get_user_inputs().user_constants.get_model_constant_rank_1_tensor("center");
};

PRISMS_PF_END_NAMESPACE