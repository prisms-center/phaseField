#ifndef dof_handler_h
#define dof_handler_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>

#include <prismspf/config.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <map>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Class that manages the deal.II DoFHandlers
 */
template <int dim>
class dofHandler
{
public:
  /**
   * \brief Constructor.
   */
  explicit dofHandler(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Destructor.
   */
  ~dofHandler();

  /**
   * \brief Initialize the DoFHandlers
   */
  void
  init(const triangulationHandler<dim>                  &triangulation_handler,
       const std::map<fieldType, dealii::FESystem<dim>> &fe_system);

  /**
   * \brief Collection of the triangulation DoFs. The number of DoFHandlers should be
   * equal to or less than the number of fields. Technically, there's a small
   * optimization we can use when multiple fields have the same constraints and
   * quadrature rule, allowing us to share the same DoFHandler. An example of this might
   * be grain growth.
   */
  std::vector<dealii::DoFHandler<dim> *> dof_handlers;

  /**
   * \brief Const copy of the dof_handlers.
   */
  std::vector<const dealii::DoFHandler<dim> *> const_dof_handlers;

private:
  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;
};

PRISMS_PF_END_NAMESPACE

#endif