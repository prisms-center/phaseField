// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/system_wide.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
DofManager<dim>::DofManager(const std::vector<FieldAttributes> &field_attributes)
{
  unsigned int num_mg_levels = 1; // Todo: upgrade to multigrid
  level_dof_handlers.resize(num_mg_levels);
  dof_handlers.resize(field_attributes.size());
  for (unsigned int field_index = 0; field_index < field_attributes.size(); ++field_index)
    {
      dof_handlers[field_index].resize(num_mg_levels, nullptr);
      for (unsigned int relative_level = 0; relative_level < num_mg_levels;
           ++relative_level)
        {
          dof_handlers[field_index][relative_level] =
            &level_dof_handlers[relative_level][field_attributes[field_index].field_type];
        }
    }
}

template <unsigned int dim>
void
DofManager<dim>::init(const TriangulationManager<dim> &triangulation_handler)
{
  for (unsigned int relative_level = 0; relative_level < level_dof_handlers.size();
       ++relative_level)
    {
      for (unsigned int rank = 0; rank < 2; ++rank)
        {
          dealii::DoFHandler<dim> &dof_handler = level_dof_handlers[relative_level][rank];
          dof_handler.reinit(triangulation_handler.get_triangulation(relative_level));
          dof_handler.distribute_dofs(SystemWide<dim, 1>::fe_systems.at(rank));
        }
    }
}

template <unsigned int dim>
void
DofManager<dim>::reinit(const TriangulationManager<dim> &triangulation_handler)
{
  init(triangulation_handler);
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

#include "core/dof_manager.inst"

PRISMS_PF_END_NAMESPACE
