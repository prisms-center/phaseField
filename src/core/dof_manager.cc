// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/system_wide.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
void
DoFManager<dim, degree>::reinit(const TriangulationManager<dim> &triangulation_manager,
                                bool                             init_mg)
{
  // reinit actual dofhandlers
  for (unsigned int rank = 0; rank < 2; ++rank)
    {
      dealii::DoFHandler<dim> &dof_handler = level_dof_handlers.at(rank);
      dof_handler.reinit(triangulation_manager.get_triangulation());
      dof_handler.distribute_dofs(SystemWide<dim, degree>::fe_systems.at(rank));
      if (init_mg)
        {
          dof_handler.distribute_mg_dofs();
        }
    }
}

template <unsigned int dim, unsigned int degree>
void
DoFManager<dim, degree>::reinit_mapping(
  const std::vector<FieldAttributes> &field_attributes)
{
  field_dof_handlers.resize(field_attributes.size(), nullptr);
  for (unsigned int field_index = 0; field_index < field_attributes.size(); ++field_index)
    {
      field_dof_handlers[field_index] = &(level_dof_handlers.at(
        static_cast<unsigned int>(field_attributes[field_index].field_type)));
    }
}

template <unsigned int dim, unsigned int degree>
const std::vector<const dealii::DoFHandler<dim> *> &
DoFManager<dim, degree>::get_field_dof_handlers() const
{
  return field_dof_handlers;
}

template <unsigned int dim, unsigned int degree>
const dealii::DoFHandler<dim> &
DoFManager<dim, degree>::get_field_dof_handler(Types::Index field_index) const
{
  return *field_dof_handlers[field_index];
}

template <unsigned int dim, unsigned int degree>
std::vector<const dealii::DoFHandler<dim> *>
DoFManager<dim, degree>::get_block_dof_handlers(
  const std::set<unsigned int> &field_indices) const
{
  std::vector<const dealii::DoFHandler<dim> *> block_dof_handlers;
  block_dof_handlers.reserve(field_indices.size());
  for (const auto &field_index : field_indices)
    {
      block_dof_handlers.push_back(field_dof_handlers[field_index]);
    }
  return block_dof_handlers;
}

template <unsigned int dim, unsigned int degree>
const std::array<dealii::DoFHandler<dim>, 2> &
DoFManager<dim, degree>::get_dof_handlers() const
{
  return level_dof_handlers;
}

template <unsigned int dim, unsigned int degree>
const dealii::DoFHandler<dim> &
DoFManager<dim, degree>::get_dof_handler(const unsigned int &rank) const
{
  return level_dof_handlers.at(rank);
}

template <unsigned int dim, unsigned int degree>
bool
DoFManager<dim, degree>::has_mg() const
{
  return level_dof_handlers[0].has_level_dofs();
}

template <unsigned int dim, unsigned int degree>
unsigned int
DoFManager<dim, degree>::num_levels() const
{
  return level_dof_handlers[0].get_triangulation().n_levels();
}

template <unsigned int dim, unsigned int degree>
dealii::types::global_dof_index
DoFManager<dim, degree>::get_total_dofs() const
{
  dealii::types::global_dof_index n_dofs = 0;
  for (const auto &dof_handler : field_dof_handlers)
    {
      n_dofs += dof_handler->n_dofs();
    }
  return n_dofs;
}

#include "core/dof_manager.inst"

PRISMS_PF_END_NAMESPACE
