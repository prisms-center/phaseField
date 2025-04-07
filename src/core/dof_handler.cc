// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <ostream>
#include <set>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
dofHandler<dim>::dofHandler(const userInputParameters<dim> &_user_inputs)
  : user_inputs(&_user_inputs)
{
  for (const auto &[index, variable] : *user_inputs->var_attributes)
    {
#ifdef ADDITIONAL_OPTIMIZATIONS
      if (user_inputs->var_attributes->at(index).duplicate_field_index !=
          numbers::invalid_index)
        {
          const_dof_handlers.push_back(
            dof_handlers.at(user_inputs->var_attributes->at(index).duplicate_field_index)
              .get());
          continue;
        }
#endif

      dof_handlers[index] = std::make_unique<dealii::DoFHandler<dim>>();
      const_dof_handlers.push_back(dof_handlers.at(index).get());
    }
}

template <int dim>
void
dofHandler<dim>::init(const triangulationHandler<dim> &triangulation_handler,
                      const std::map<fieldType, dealii::FESystem<dim>> &fe_system)
{
  // TODO (landinjm): Include multigrid degrees of freedom.
  unsigned int n_dofs = 0;
  for (const auto &[index, variable] : *user_inputs->var_attributes)
    {
#ifdef ADDITIONAL_OPTIMIZATIONS
      if (user_inputs->var_attributes->at(index).duplicate_field_index !=
          numbers::invalid_index)
        {
          n_dofs +=
            dof_handlers.at(user_inputs->var_attributes->at(index).duplicate_field_index)
              ->n_dofs();
          continue;
        }
#endif
      dof_handlers.at(index)->reinit(triangulation_handler.get_triangulation());
      dof_handlers.at(index)->distribute_dofs(fe_system.at(variable.field_type));

      n_dofs += dof_handlers.at(index)->n_dofs();
    }

  // Check if multigrid is enabled
  std::set<unsigned int> fields_with_multigrid;
  for (const auto &[index, linear_solver_parameters] :
       user_inputs->linear_solve_parameters.linear_solve)
    {
      if (linear_solver_parameters.preconditioner == preconditionerType::GMG)
        {
          has_multigrid = true;
          fields_with_multigrid.insert(index);
        }
    }

  // If we don't have multigrid print relevant info and return early
  if (!has_multigrid)
    {
      // TODO (landinjm): Print other useful information in debug mode
      conditionalOStreams::pout_base()
        << "  number of degrees of freedom: " << n_dofs << "\n"
        << std::flush;
      return;
    }

  // TODO (landinjm): Add optimizations for shared DoFHandlers

  const unsigned int min_level      = triangulation_handler.get_mg_min_level();
  const unsigned int max_level      = triangulation_handler.get_mg_max_level();
  unsigned int       n_dofs_with_mg = n_dofs;
  for (const auto &[index, variable] : *user_inputs->var_attributes)
    {
      if (fields_with_multigrid.find(index) == fields_with_multigrid.end())
        {
          continue;
        }

      mg_dof_handlers.emplace(index,
                              dealii::MGLevelObject<dealii::DoFHandler<dim>>(min_level,
                                                                             max_level));

      for (unsigned int level = min_level; level <= max_level; ++level)
        {
          mg_dof_handlers.at(index)[level].reinit(
            triangulation_handler.get_mg_triangulation(level));
          mg_dof_handlers.at(index)[level].distribute_dofs(
            fe_system.at(variable.field_type));
          n_dofs_with_mg += mg_dof_handlers.at(index)[level].n_dofs();
        }
    }

  // TODO (landinjm): Print other useful information in debug mode
  conditionalOStreams::pout_base()
    << "  number of degrees of freedom: " << n_dofs << "\n"
    << "    with multigrid levels: " << n_dofs_with_mg << "\n"
    << std::flush;
}

template <int dim>
const std::vector<const dealii::DoFHandler<dim> *> &
dofHandler<dim>::get_dof_handlers() const
{
  Assert(const_dof_handlers.size() == user_inputs->var_attributes->size(),
         dealii::ExcNotInitialized());
  return const_dof_handlers;
}

template <int dim>
const std::map<unsigned int, dealii::MGLevelObject<dealii::DoFHandler<dim>>> &
dofHandler<dim>::get_mg_dof_handlers() const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());
  Assert(!mg_dof_handlers.empty(),
         dealii::ExcMessage("The multigrid dof handler map is empty."));
  return mg_dof_handlers;
}

template <int dim>
const std::vector<const dealii::DoFHandler<dim> *> &
dofHandler<dim>::get_mg_dof_handlers(unsigned int level) const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());
  Assert(!mg_dof_handlers.empty(),
         dealii::ExcMessage("The multigrid dof handler map is empty."));
  Assert(level >= const_mg_dof_handlers.min_level() &&
           level <= const_mg_dof_handlers.max_level(),
         dealii::ExcIndexRange(level,
                               const_mg_dof_handlers.min_level(),
                               const_mg_dof_handlers.max_level()));
  return const_mg_dof_handlers[level];
}

template <int dim>
const dealii::DoFHandler<dim> &
dofHandler<dim>::get_mg_dof_handler(unsigned int index, unsigned int level) const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());
  Assert(!mg_dof_handlers.empty(),
         dealii::ExcMessage("The multigrid dof handler map is empty."));
  Assert(mg_dof_handlers.find(index) != mg_dof_handlers.end(),
         dealii::ExcMessage("The multigrid dof handler map does not contain the index."));
  Assert(level >= mg_dof_handlers.at(index).min_level() &&
           level <= mg_dof_handlers.at(index).max_level(),
         dealii::ExcIndexRange(level,
                               mg_dof_handlers.at(index).min_level(),
                               mg_dof_handlers.at(index).max_level()));
  return mg_dof_handlers.at(index)[level];
}

INSTANTIATE_UNI_TEMPLATE(dofHandler)

PRISMS_PF_END_NAMESPACE
