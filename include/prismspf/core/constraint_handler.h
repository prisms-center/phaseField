#ifndef constraint_handler_h
#define constraint_handler_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/affine_constraints.h>

#include <prismspf/config.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief The class handles the generation and application of boundary conditions based on
 * the user-inputs.
 */
template <int dim>
class constraintHandler
{
public:
  /**
   * \brief Constructor.
   */
  explicit constraintHandler(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Getter function for the constraints.
   */
  [[nodiscard]] std::vector<const dealii::AffineConstraints<double> *>
  get_constraints();

  /**
   * \brief Getter function for the constraint of an index (constant reference).
   */
  [[nodiscard]] const dealii::AffineConstraints<double> &
  get_constraint(const unsigned int &index) const;

  /**
   * \brief Make constraints based on the inputs of the construct.
   */
  void
  make_constraints(const dealii::Mapping<dim>                   &mapping,
                   const std::vector<dealii::DoFHandler<dim> *> &dof_handlers);

private:
  /**
   * \brief Make the constrainst for a single index.
   */
  void
  make_constraint(const dealii::Mapping<dim>    &mapping,
                  const dealii::DoFHandler<dim> &dof_handler,
                  const unsigned int            &index);

  /**
   * \brief Set the dirichlet constraint for the pinned point.
   */
  void
  set_pinned_point(const dealii::DoFHandler<dim> &dof_handler, const unsigned int &index);

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief Constraints.
   */
  std::vector<dealii::AffineConstraints<double>> constraints;
};

PRISMS_PF_END_NAMESPACE

#endif