#ifndef matrix_free_handler_h
#define matrix_free_handler_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <core/user_inputs/user_input_parameters.h>
#include <vector>

/**
 * \brief This class handlers the management and access of the matrix-free objects.
 */
template <int dim, typename number = double>
class matrixfreeHandler
{
public:
  /**
   * \brief Constructor.
   */
  matrixfreeHandler(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Reinitialize the matrix-free object with the same quad rule.
   */
  void
  reinit(const dealii::Mapping<dim>                                   &mapping,
         const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
         const std::vector<const dealii::AffineConstraints<number> *> &constraint,
         const Quadrature<1>                                          &quad);

  /**
   * \brief Reinitialize the matrix-free object with the different quad rule.
   */
  void
  reinit(const dealii::Mapping<dim>                                   &mapping,
         const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
         const std::vector<const dealii::AffineConstraints<number> *> &constraint,
         const std::vector<Quadrature<1>>                             &quad);

  /**
   * \brief Getter function for the matrix-free object (shared ptr).
   */
  std::shared_ptr<dealii::MatrixFree<dim, number>>
  get_matrix_free();

private:
  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief Matrix-free object that collects data to be used in cell loop operations.
   */
  std::shared_ptr<dealii::MatrixFree<dim, number>> matrix_free_object;

  /**
   * \brief Additional data scheme
   */
  typename dealii::MatrixFree<dim, number>::AdditionalData additional_data;
};

template <int dim, typename number>
matrixfreeHandler<dim, number>::matrixfreeHandler(
  const userInputParameters<dim> &_user_inputs)
  : user_inputs(_user_inputs)
  , matrix_free_object(std::make_shared<dealii::MatrixFree<dim, number>>())
{
  additional_data.tasks_parallel_scheme =
    dealii::MatrixFree<dim, number>::AdditionalData::none;

  // This should be done according to the residual flags to prevent excess data being
  // evaluated for the shape functions. Note that this applies to all PDEs.
  additional_data.mapping_update_flags =
    (update_values | update_gradients | update_JxW_values | update_quadrature_points);
}

template <int dim, typename number>
inline void
matrixfreeHandler<dim, number>::reinit(
  const dealii::Mapping<dim>                                   &mapping,
  const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
  const std::vector<const dealii::AffineConstraints<number> *> &constraint,
  const Quadrature<1>                                          &quad)
{
  matrix_free_object->reinit(mapping, dof_handler, constraint, quad, additional_data);
}

template <int dim, typename number>
inline void
matrixfreeHandler<dim, number>::reinit(
  const dealii::Mapping<dim>                                   &mapping,
  const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
  const std::vector<const dealii::AffineConstraints<number> *> &constraint,
  const std::vector<Quadrature<1>>                             &quad)
{
  matrix_free_object->reinit(mapping, dof_handler, constraint, quad, additional_data);
}

template <int dim, typename number>
inline std::shared_ptr<dealii::MatrixFree<dim, number>>
matrixfreeHandler<dim, number>::get_matrix_free()
{
  return matrix_free_object;
}

#endif