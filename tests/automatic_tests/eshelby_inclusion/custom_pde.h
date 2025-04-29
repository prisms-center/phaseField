// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

//#include <core/matrixFreePDE.h> //JM moved; maybe delete

//JM copied includes from precipitate app
#pragma once

#include <prismspf/core/pde_operator.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>


// using namespace dealii; //JM gone in 3.0
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

template <unsigned int dim, unisgned int degree, typename number> //JM unsigned ints and new typename
//class customPDE : public MatrixFreePDE<dim, degree> //JM new class
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

  //JM constructor changes
  /*
  // Constructor
  customPDE(userInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
    , userInputs(_userInputs) {};
  */

  /**
  * \brief Constructor.
  */
  explicit customPDE(const userInputParameters<dim> &_user_inputs)
  : PDEOperator<dim, degree, number>(_user_inputs)
  {}
  //JM leaving this for now, sets IC's + BC's, prob changed in 3.0
  // Function to set the initial conditions (in ICs_and_BCs.h)
  void
  setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                      [[maybe_unused]] const unsigned int index,
                      [[maybe_unused]] double            &scalar_IC,
                      [[maybe_unused]] Vector<double>    &vector_IC) override;

  // Function to set the non-uniform Dirichlet boundary conditions (in
  // ICs_and_BCs.h)
  void
  setNonUniformDirichletBCs([[maybe_unused]] const Point<dim>  &p,
                            [[maybe_unused]] const unsigned int index,
                            [[maybe_unused]] const unsigned int direction,
                            [[maybe_unused]] const double       time,
                            [[maybe_unused]] double            &scalar_BC,
                            [[maybe_unused]] Vector<double>    &vector_BC) override;

private:
#include <core/typeDefs.h>
  /**
   * \brief User-implemented class for the RHS of explicit equations.
   */
  //const userInputParameters<dim> userInputs;

  // Function to set the RHS of the governing equations for explicit time
  // dependent equations (in equations.h)
  //JM old explicit RHS
  /*
  void
  explicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;
  */
  //JM new explicit RHS
  void
  compute_explicit_RHS(variableContainer<dim, degree, number, dealii::VectorizedArray<double>> 
                                                                    &variable_list,
                       const dealii::Point<dim, dealii::VectorizedArray<number>>
                         &q_point_loc
                       const dealii::VectorizedArray<number> element_volume) const override;

  //JM old nonexplicit RHS
  // Function to set the RHS of the governing equations for all other equations
  // (in equations.h)
  /*
  void
  nonExplicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;
  */
  /**
    * \brief User-implemented class for the RHS of nonexplicit equations.
    */
  //JM new nonexplicit RHS
  void
  compute_nonexplicit_RHS(
    variableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
    types::index current_index = numbers::invalid_index) const override;
 
    //JM note, eshelby has no existing differentiation between explicit and nonexplicit LHS
    // Function to set the LHS of the governing equations (in equations.h)
  void
  equationLHS(
    variableContainer<dim, degree, VectorizedArray<number>>
                                                              &variable_list,
    const Point<dim, VectorizedArray<double>> q_point_loc,
    const VectorizedArray<double> element_volume) const override;

//JM Slated to overhaul, for now it is commented out
// Function to set postprocessing expressions (in postprocess.h)
/*
#ifdef POSTPROCESS_FILE_EXISTS
  void
  postProcessedFields(
    [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
      &variable_list,
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &pp_variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;
#endif
*/
//JM I don't think this exists anymore
/*
// Function to set the nucleation probability (in nucleation.h)
#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability([[maybe_unused]] variableValueContainer variable_value,
                           [[maybe_unused]] double                 dV) const override;
#endif
*/

  // ================================================================
  // Methods specific to this subclass
  // ================================================================

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  constexpr static unsigned int  CIJ_tensor_size = 2 * dim - 1 + dim / 3; //JM modified to constexpr
  //JM syntax update for stiffness tensor and other constants
  dealii::Tensor<2, CIJ_tensor_size, number> CIJ = 
    this->get_user_inputs().user_constants.get_model_constant_elasticity_tensor("CIJ");

  number incRadius = this->get_user_inputs().user_constants.get_model_constant_double("incRadius");
  number poisson = this->get_user_inputs().user_constants.get_model_constant_double("poisson");
  dealii::Tensor<1, dim, number> center =
  this->get_user_inputs().user_constants.get_model_constant_rank_1_tensor("center");

  // ================================================================
};