// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <memory>
#include <ostream>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
DofHandler<dim>::DofHandler(const UserInputParameters<dim> &_user_inputs,
                            const MGInfo<dim>              &mg_info)
  : user_inputs(&_user_inputs)
{
  for (const auto &[index, variable] : user_inputs->get_variable_attributes())
    {
#ifdef ADDITIONAL_OPTIMIZATIONS
      // TODO (landinjm): This relies on the fact the entry has already been created. Add
      // an assertion
      if (user_inputs->get_variable_attributes().at(index).get_duplicate_field_index() !=
          Numbers::invalid_index)
        {
          const_dof_handlers.push_back(dof_handlers
                                         .at(user_inputs->get_variable_attributes()
                                               .at(index)
                                               .get_duplicate_field_index())
                                         .get());
          continue;
        }
#endif

      dof_handlers[index] = std::make_unique<dealii::DoFHandler<dim>>();
      const_dof_handlers.push_back(dof_handlers.at(index).get());
    }
  // If we don't have multigrid, we can return early
  if (!mg_info.has_multigrid())
    {
      return;
    }
  has_multigrid = true;

  // Go through all fields that have multigrid levels and create the DoFHandlers
  global_min_level                    = mg_info.get_mg_min_level();
  const unsigned int global_max_level = mg_info.get_mg_max_level();
  const_mg_dof_handlers.resize(global_max_level - global_min_level + 1);
  for (const auto &[index, dependency, min_level] : mg_info.get_lhs_fields())
    {
      const unsigned int relative_level = min_level - global_min_level;
#ifdef ADDITIONAL_OPTIMIZATIONS
      // TODO (landinjm): This relies on the fact the entry has already been created. Add
      // an assertion

      // TODO (landinjm): Duplicate field index is not supported yet because we have
      // to handle the case where the duplicate field has different min multigrid
      // levels.
      AssertThrow(false,
                  FeatureNotImplemented("Additional optimizations for multigrid fields"));
      // TODO (landinjm): Add n_mg_dofs

#endif
      if (mg_dof_handlers.contains(index))
        {
          // TODO (landinjm): Small edge case where the Change and Normal term have
          // different min levels.
          Assert(mg_dof_handlers.at(index).min_level() == min_level,
                 dealii::ExcMessage("The minimum multigrid level for index " +
                                    std::to_string(index) +
                                    " is not the same as the one in the MGInfo class."));
          continue;
        }

      mg_dof_handlers[index] =
        dealii::MGLevelObject<std::unique_ptr<dealii::DoFHandler<dim>>>(min_level,
                                                                        global_max_level);
      for (unsigned int level = relative_level; level < const_mg_dof_handlers.size();
           level++)
        {
          mg_dof_handlers[index][level + global_min_level] =
            std::make_unique<dealii::DoFHandler<dim>>();
          const_mg_dof_handlers[level].push_back(
            mg_dof_handlers[index][level + global_min_level].get());
        }
    }
}

template <unsigned int dim>
void
DofHandler<dim>::init(const TriangulationHandler<dim> &triangulation_handler,
                      const std::map<FieldType, dealii::FESystem<dim>> &fe_system,
                      const MGInfo<dim>                                &mg_info)
{
  unsigned int n_dofs = 0;
  for (const auto &[index, variable] : user_inputs->get_variable_attributes())
    {
#ifdef ADDITIONAL_OPTIMIZATIONS
      if (user_inputs->get_variable_attributes().at(index).get_duplicate_field_index() !=
          Numbers::invalid_index)
        {
          n_dofs += dof_handlers
                      .at(user_inputs->get_variable_attributes()
                            .at(index)
                            .get_duplicate_field_index())
                      ->n_dofs();
          continue;
        }
#endif
      dof_handlers.at(index)->reinit(triangulation_handler.get_triangulation());
      dof_handlers.at(index)->distribute_dofs(fe_system.at(variable.get_field_type()));

      n_dofs += dof_handlers.at(index)->n_dofs();
    }

  // If we don't have multigrid print relevant info and return early
  if (!has_multigrid)
    {
      // TODO (landinjm): Print other useful information in debug mode
      ConditionalOStreams::pout_base()
        << "  number of degrees of freedom: " << n_dofs << "\n"
        << std::flush;
      return;
    }

  unsigned int n_dofs_with_mg = n_dofs;

  // Go through all fields that have multigrid levels and reinit the DoFHandlers
  for (const auto &[index, dependency, min_level] : mg_info.get_lhs_fields())
    {
#ifdef ADDITIONAL_OPTIMIZATIONS
      AssertThrow(false,
                  FeatureNotImplemented("Additional optimizations for multigrid fields"));
#endif
      for (unsigned int level = min_level; level <= mg_info.get_mg_max_level(); ++level)
        {
          mg_dof_handlers.at(index)[level]->reinit(
            triangulation_handler.get_mg_triangulation(level));
          mg_dof_handlers.at(index)[level]->distribute_dofs(fe_system.at(
            user_inputs->get_variable_attributes().at(index).get_field_type()));
          n_dofs_with_mg += mg_dof_handlers.at(index)[level]->n_dofs();
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
DofHandler<dim>::reinit(const TriangulationHandler<dim> &triangulation_handler,
                        const std::map<FieldType, dealii::FESystem<dim>> &fe_system,
                        const MGInfo<dim>                                &mg_info)
{
  this->init(triangulation_handler, fe_system, mg_info);
}

template <unsigned int dim>
const std::vector<const dealii::DoFHandler<dim> *> &
DofHandler<dim>::get_dof_handlers() const
{
  Assert(const_dof_handlers.size() == user_inputs->get_variable_attributes().size(),
         dealii::ExcNotInitialized());
  return const_dof_handlers;
}

template <unsigned int dim>
const std::vector<const dealii::DoFHandler<dim> *> &
DofHandler<dim>::get_mg_dof_handlers(unsigned int level) const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());
  Assert(!const_mg_dof_handlers.empty(), dealii::ExcNotInitialized());
  Assert((level - global_min_level) < const_mg_dof_handlers.size(),
         dealii::ExcMessage("The requested level is out of range."));
  return const_mg_dof_handlers[level - global_min_level];
}

INSTANTIATE_UNI_TEMPLATE(DofHandler)

PRISMS_PF_END_NAMESPACE
