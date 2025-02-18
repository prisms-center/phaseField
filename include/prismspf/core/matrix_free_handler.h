#ifndef matrix_free_handler_h
#define matrix_free_handler_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/config.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <memory>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

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
  reinit(const dealii::Mapping<dim>              &mapping,
         const dealii::DoFHandler<dim>           &dof_handler,
         const dealii::AffineConstraints<number> &constraint,
         const dealii::Quadrature<1>             &quad);

  /**
   * \brief Reinitialize the matrix-free object with the same quad rule.
   */
  void
  reinit(const dealii::Mapping<dim>                                   &mapping,
         const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
         const std::vector<const dealii::AffineConstraints<number> *> &constraint,
         const dealii::Quadrature<1>                                  &quad);

  /**
   * \brief Reinitialize the matrix-free object with the different quad rule.
   */
  void
  reinit(const dealii::Mapping<dim>                                   &mapping,
         const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
         const std::vector<const dealii::AffineConstraints<number> *> &constraint,
         const std::vector<dealii::Quadrature<1>>                     &quad);

  /**
   * \brief Getter function for the matrix-free object (shared ptr).
   */
  [[nodiscard]] std::shared_ptr<dealii::MatrixFree<dim, number>>
  get_matrix_free() const;

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

PRISMS_PF_END_NAMESPACE

#endif