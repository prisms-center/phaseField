// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/system_wide.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
DoFManager<dim, degree>::DoFManager(
  const std::vector<FieldAttributes> &field_attributes,
  const TriangulationManager<dim>    &triangulation_manager)
{
  init(1);
  reinit(triangulation_manager);
  reinit_mapping(field_attributes);
}

template <unsigned int dim, unsigned int degree>
void
DoFManager<dim, degree>::init(unsigned int num_levels)
{
  level_dof_handlers = std::vector<std::array<dealii::DoFHandler<dim>, 2>>(num_levels);
  field_dof_handlers.clear();
}

template <unsigned int dim, unsigned int degree>
void
DoFManager<dim, degree>::reinit(const TriangulationManager<dim> &triangulation_manager)
{
  for (unsigned int relative_level = 0; relative_level < level_dof_handlers.size();
       ++relative_level)
    {
      // reinit actual dofhandlers
      for (unsigned int rank = 0; rank < 2; ++rank)
        {
          dealii::DoFHandler<dim> &dof_handler =
            level_dof_handlers[relative_level].at(rank);
          dof_handler.reinit(triangulation_manager.get_triangulation(relative_level));
          dof_handler.distribute_dofs(SystemWide<dim, degree>::fe_systems.at(rank));
        }
    }
}

template <unsigned int dim, unsigned int degree>
void
DoFManager<dim, degree>::reinit_mapping(
  const std::vector<FieldAttributes> &field_attributes)
{
  const unsigned int num_levels = level_dof_handlers.size();
  field_dof_handlers.resize(num_levels);
  for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
    {
      field_dof_handlers[relative_level].resize(field_attributes.size(), nullptr);
      for (unsigned int field_index = 0; field_index < field_attributes.size();
           ++field_index)
        {
          field_dof_handlers[relative_level][field_index] =
            &(level_dof_handlers[relative_level].at(
              static_cast<unsigned int>(field_attributes[field_index].field_type)));
        }
    }
}

template <unsigned int dim, unsigned int degree>
const std::vector<std::vector<const dealii::DoFHandler<dim> *>> &
DoFManager<dim, degree>::get_field_dof_handlers_levels() const
{
  return field_dof_handlers;
}

template <unsigned int dim, unsigned int degree>
const std::vector<const dealii::DoFHandler<dim> *> &
DoFManager<dim, degree>::get_field_dof_handlers(unsigned int relative_level) const
{
  return field_dof_handlers[relative_level];
}

template <unsigned int dim, unsigned int degree>
const dealii::DoFHandler<dim> &
DoFManager<dim, degree>::get_field_dof_handler(Types::Index field_index,
                                               unsigned int relative_level) const
{
  return *field_dof_handlers[relative_level][field_index];
}

template <unsigned int dim, unsigned int degree>
const std::vector<std::array<dealii::DoFHandler<dim>, 2>> &
DoFManager<dim, degree>::get_dof_handlers_levels() const
{
  return level_dof_handlers;
}

template <unsigned int dim, unsigned int degree>
const std::array<dealii::DoFHandler<dim>, 2> &
DoFManager<dim, degree>::get_dof_handlers(unsigned int relative_level) const
{
  return level_dof_handlers[relative_level];
}

template <unsigned int dim, unsigned int degree>
const dealii::DoFHandler<dim> &
DoFManager<dim, degree>::get_dof_handler(const unsigned int &rank,
                                         unsigned int        relative_level) const
{
  return level_dof_handlers[relative_level].at(rank);
}

template <unsigned int dim, unsigned int degree>
dealii::types::global_dof_index
DoFManager<dim, degree>::get_total_dofs() const
{
  dealii::types::global_dof_index n_dofs = 0;
  for (const auto &dof_handler : field_dof_handlers[0])
    {
      n_dofs += dof_handler->n_dofs();
    }
  return n_dofs;
}

#include "core/dof_manager.inst"

PRISMS_PF_END_NAMESPACE
