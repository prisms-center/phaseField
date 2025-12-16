// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/exceptions.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
DofManager<dim>::DofManager(const std::vector<FieldAttributes> &field_attributes,
                            const std::set<SolveGroup>         &solve_groups)
{
  dof_handlers.resize(field_attributes.size());
  for (const auto &solve_group : solve_groups)
    {
      /* TODO (fractalsbyx) figure out mg depth. Careful. Aux fields inherit this from
       * primary fields, as well as any dependencies*/
      unsigned int num_mg_levels = 1;
      for (const auto &field_index : solve_group.field_indices)
        {
          dof_handlers[field_index].resize(num_mg_levels, nullptr);
          for (unsigned int relative_level = 0; relative_level < num_mg_levels;
               relative_level++)
            {
              dof_handlers[field_index][relative_level] =
                std::make_shared<dealii::DoFHandler<dim>>();
            }
        }
    }
}

template <unsigned int dim>
void
DofManager<dim>::init(const TriangulationManager<dim>    &triangulation_handler,
                      const std::vector<FieldAttributes> &field_attributes)
{
  unsigned int n_dofs         = 0;
  unsigned int n_dofs_with_mg = 0;
  for (unsigned int field_index = 0; field_index < field_attributes.size(); ++field_index)
    {
      for (unsigned int relative_level = 0;
           relative_level < dof_handlers[field_index].size();
           ++relative_level)
        {
          std::shared_ptr<dealii::DoFHandler<dim>> &dof_handler =
            dof_handlers[field_index][relative_level];
          dof_handler->reinit(triangulation_handler.get_triangulation(relative_level));
          dof_handler->distribute_dofs(
            fe_systems.at(field_attributes[field_index].field_type));
          n_dofs_with_mg += dof_handler->n_dofs();
          n_dofs += bool(relative_level) ? 0 : dof_handler->n_dofs();
        }
    }

  // TODO (landinjm): Print other useful information in debug mode
  ConditionalOStreams::pout_base()
    << "  number of degrees of freedom: " << n_dofs << "\n"
    << "    with multigrid levels: " << n_dofs_with_mg << "\n"
    << std::flush;
}

template <unsigned int dim>
void
DofManager<dim>::reinit(const TriangulationManager<dim>    &triangulation_handler,
                        const std::vector<FieldAttributes> &field_attributes)
{
  init(triangulation_handler, field_attributes);
}

template <unsigned int dim>
const std::vector<const dealii::DoFHandler<dim> *> &
DofManager<dim>::get_dof_handlers(const std::set<unsigned int> &field_indices,
                                  unsigned int                  relative_level) const
{
  std::vector<const dealii::DoFHandler<dim> *> selected_dof_handlers;
  selected_dof_handlers.reserve(field_indices.size());
  for (const auto index : field_indices)
    {
      selected_dof_handlers.push_back(dof_handlers[index][relative_level].get());
    }
  return selected_dof_handlers;
}

// #include "core/dof_manager.inst"

PRISMS_PF_END_NAMESPACE
