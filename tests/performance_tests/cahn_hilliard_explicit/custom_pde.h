// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <cmath>

PRISMS_PF_BEGIN_NAMESPACE

const unsigned int n_copies = 1;

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
};

inline void
CustomAttributeLoader::load_variable_attributes()
{
  for (unsigned int i = 0; i < n_copies; i++)
    {
      std::string field_name     = "phi" + std::to_string(i);
      std::string aux_field_name = "mu" + std::to_string(i);

      set_variable_name(i, field_name);
      set_variable_type(i, FieldInfo::TensorRank::Scalar);
      set_variable_equation_type(i, ExplicitTimeDependent);

      set_dependencies_value_term_rhs(i, field_name);
      set_dependencies_gradient_term_rhs(i, "grad(" + aux_field_name + ")");

      set_variable_name(n_copies + i, aux_field_name);
      set_variable_type(n_copies + i, FieldInfo::TensorRank::Scalar);
      set_variable_equation_type(n_copies + i, Auxiliary);

      set_dependencies_value_term_rhs(n_copies + i, field_name);
      set_dependencies_gradient_term_rhs(n_copies + i, "grad(" + field_name + ")");
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::set_initial_condition(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{
  if (index < n_copies)
    {
      double center[12][3] = {
        {0.1, 0.3,  0},
        {0.8, 0.7,  0},
        {0.5, 0.2,  0},
        {0.4, 0.4,  0},
        {0.3, 0.9,  0},
        {0.8, 0.1,  0},
        {0.9, 0.5,  0},
        {0.0, 0.1,  0},
        {0.1, 0.6,  0},
        {0.5, 0.6,  0},
        {1,   1,    0},
        {0.7, 0.95, 0}
      };
      double rad[12] = {12, 14, 19, 16, 11, 12, 17, 15, 20, 10, 11, 14};
      double dist    = 0.0;
      for (unsigned int i = 0; i < 12; i++)
        {
          dist = 0.0;
          for (unsigned int dir = 0; dir < dim; dir++)
            {
              dist +=
                (point[dir] - center[i][dir] * 100) * (point[dir] - center[i][dir] * 100);
            }
          dist = std::sqrt(dist);

          scalar_value += 0.5 * (1.0 - std::tanh((dist - rad[i]) / 1.5));
        }
      scalar_value = std::min(scalar_value, static_cast<number>(1.0));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::set_nonuniform_dirichlet(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &boundary_id,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  const dealii::VectorizedArray<number> &element_volume,
  Types::Index                           solve_block) const
{
  for (unsigned int i = 0; i < n_copies; i++)
    {
      ScalarValue n   = variable_list.template get_value<ScalarValue>(i);
      ScalarGrad  mux = variable_list.template get_gradient<ScalarGrad>(n_copies + i);

      ScalarValue eq_n  = n;
      ScalarGrad  eqx_n = -get_timestep() * mux;

      variable_list.set_value_term(i, eq_n);
      variable_list.set_gradient_term(i, eqx_n);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  const dealii::VectorizedArray<number> &element_volume,
  Types::Index                           solve_block,
  [[maybe_unused]] Types::Index          index) const
{
  for (unsigned int i = 0; i < n_copies; i++)
    {
      if (index == n_copies + i)
        {
          ScalarValue n  = variable_list.template get_value<ScalarValue>(i);
          ScalarGrad  nx = variable_list.template get_gradient<ScalarGrad>(i);

          ScalarValue fnV    = 4.0 * n * (n - 1.0) * (n - 0.5);
          ScalarValue eq_mu  = fnV;
          ScalarGrad  eqx_mu = 2.0 * nx;

          variable_list.set_value_term(n_copies + i, eq_mu);
          variable_list.set_gradient_term(n_copies + i, eqx_mu);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  const dealii::VectorizedArray<number> &element_volume,
  Types::Index                           solve_block,
  [[maybe_unused]] Types::Index          index) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  const dealii::VectorizedArray<number> &element_volume,
  Types::Index                           solve_block) const
{}

INSTANTIATE_TRI_TEMPLATE(CustomPDE)

PRISMS_PF_END_NAMESPACE