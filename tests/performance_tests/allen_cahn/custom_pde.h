// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/nonuniform_dirichlet.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <algorithm>
#include <cmath>

PRISMS_PF_BEGIN_NAMESPACE

const unsigned int n_copies = 64;

/**
 * \brief This is a derived class of `MatrixFreeOperator` where the user implements their
 * PDEs.
 *
 * \tparam dim The number of dimensions in the problem.
 * \tparam degree The polynomial degree of the shape functions.
 * \tparam number Datatype to use. Either double or float.
 */
template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public MatrixFreeOperator<dim, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using ScalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using VectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;

  /**
   * \brief Constructor for concurrent solves.
   */
  CustomPDE(const UserInputParameters<dim>                   &_user_inputs,
            const std::map<unsigned int, VariableAttributes> &subset_attributes)
    : MatrixFreeOperator<dim, degree, number>(_user_inputs, subset_attributes)
  {}

  /**
   * \brief Constructor for single solves.
   */
  CustomPDE(const UserInputParameters<dim>                   &_user_inputs,
            const unsigned int                               &_current_index,
            const std::map<unsigned int, VariableAttributes> &subset_attributes)
    : MatrixFreeOperator<dim, degree, number>(_user_inputs,
                                              _current_index,
                                              subset_attributes)
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
  compute_nonexplicit_rhs(VariableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, dealii::VectorizedArray<number>>
                            &q_point_loc) const override;

  /**
   * \brief User-implemented class for the LHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_lhs(VariableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, dealii::VectorizedArray<number>>
                            &q_point_loc) const override;

  /**
   * \brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  void
  compute_postprocess_explicit_rhs(
    VariableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
    const override;
};

inline void
CustomAttributeLoader::load_variable_attributes()
{
  for (unsigned int i = 0; i < n_copies; i++)
    {
      std::string field_name = "phi" + std::to_string(i);

      set_variable_name(i, field_name);
      set_variable_type(i, Scalar);
      set_variable_equation_type(i, ExplicitTimeDependent);

      set_dependencies_value_term_rhs(i, field_name);
      set_dependencies_gradient_term_rhs(i, "grad(" + field_name + ")");
    }
}

template <unsigned int dim>
inline void
customInitialCondition<dim>::set_initial_condition(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] double                   &scalar_value,
  [[maybe_unused]] double                   &vector_component_value) const
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
  scalar_value = std::min(scalar_value, 1.0);
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
inline void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  for (unsigned int i = 0; i < n_copies; i++)
    {
      ScalarValue n  = variable_list.get_scalar_value(i);
      ScalarGrad  nx = variable_list.get_scalar_gradient(i);

      ScalarValue fnV   = (4.0 * n * (n - 1.0) * (n - 0.5));
      ScalarValue eq_n  = (n - (this->get_timestep() * fnV));
      ScalarGrad  eqx_n = (-this->get_timestep() * 2.0 * nx);

      variable_list.set_scalar_value_term(i, eq_n);
      variable_list.set_scalar_gradient_term(i, eqx_n);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

template <unsigned int dim, unsigned int degree, typename number>
inline void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

template <unsigned int dim, unsigned int degree, typename number>
inline void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

INSTANTIATE_TRI_TEMPLATE(CustomPDE)

PRISMS_PF_END_NAMESPACE