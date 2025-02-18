#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <map>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
dofHandler<dim>::dofHandler(const userInputParameters<dim> &_user_inputs)
  : user_inputs(_user_inputs)
{
  dof_handlers.resize(user_inputs.var_attributes.size());
  for (const auto &[index, variable] : user_inputs.var_attributes)
    {
      dof_handlers.at(index) = new dealii::DoFHandler<dim>();
    }
  for (auto &dof_handler : dof_handlers)
    {
      const_dof_handlers.push_back(dof_handler);
    }
}

template <int dim>
dofHandler<dim>::~dofHandler()
{
  for (auto dof_handler : dof_handlers)
    {
      delete dof_handler;
    }
  dof_handlers.clear();
}

template <int dim>
void
dofHandler<dim>::init(const triangulationHandler<dim> &triangulation_handler,
                      const std::map<fieldType, dealii::FESystem<dim>> &fe_system)
{
  // TODO: Fix so we can make unique instances of dof handlers.
  unsigned int n_dofs = 0;
  for (const auto &[index, variable] : user_inputs.var_attributes)
    {
      dof_handlers.at(index)->reinit(triangulation_handler.get_triangulation());
      dof_handlers.at(index)->distribute_dofs(fe_system.at(variable.field_type));

      n_dofs += dof_handlers.at(index)->n_dofs();
    }
  conditionalOStreams::pout_base() << "number of degrees of freedom: " << n_dofs << "\n";
  conditionalOStreams::pout_summary()
    << "  number of degrees of freedom: " << n_dofs << "\n"
    << std::flush;
}

INSTANTIATE_UNI_TEMPLATE(dofHandler)

PRISMS_PF_END_NAMESPACE