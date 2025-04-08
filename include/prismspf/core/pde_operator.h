#pragma once

#include <deal.II/base/point.h>

#include <prismspf/core/variable_container.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class contains the user implementation of each PDE operator.
 */
template <int dim, int degree, typename number>
class PDEOperator
{
public:
  using size_type = dealii::VectorizedArray<number>;

  /**
   * \brief Constructor.
   */
  explicit PDEOperator(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Destructor.
   */
  virtual ~PDEOperator() = default;

  /**
   * \brief User-implemented class for the RHS of explicit equations.
   */
  virtual void
  compute_explicit_RHS(variableContainer<dim, degree, number> &variable_list,
                       const dealii::Point<dim, size_type>    &q_point_loc) const = 0;

  /**
   * \brief User-implemented class for the RHS of nonexplicit equations.
   */
  virtual void
  compute_nonexplicit_RHS(variableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, size_type>    &q_point_loc,
                          types::index current_index = numbers::invalid_index) const = 0;

  /**
   * \brief User-implemented class for the LHS of nonexplicit equations.
   */
  virtual void
  compute_nonexplicit_LHS(variableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, size_type>    &q_point_loc,
                          types::index current_index = numbers::invalid_index) const = 0;

  /**
   * \brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  virtual void
  compute_postprocess_explicit_RHS(
    variableContainer<dim, degree, number> &variable_list,
    const dealii::Point<dim, size_type>    &q_point_loc) const = 0;

  /**
   * \brief Get the user inputs (constant reference).
   */
  [[nodiscard]] const userInputParameters<dim> &
  get_user_inputs() const;

  /**
   * \brief Get the timestep (copy).
   */
  [[nodiscard]] number
  get_timestep() const;

private:
  /**
   * \brief The user-inputs.
   */
  const userInputParameters<dim> *user_inputs;
};

PRISMS_PF_END_NAMESPACE