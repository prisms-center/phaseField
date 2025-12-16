// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/triangulation_manager.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <map>
#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Class that manages the deal.II DoFHandlers
 */
template <unsigned int dim>
class DofManager
{
public:
  /**
   * @brief Constructor.
   */
  DofManager(const std::vector<FieldAttributes> &field_attributes,
             const std::set<SolveGroup>         &solve_groups);

  /**
   * @brief Initialize the DoFHandlers
   */
  void
  init(const TriangulationManager<dim>    &triangulation_handler,
       const std::vector<FieldAttributes> &field_attributes);

  /**
   * @brief Reinitialize the DoFHandlers
   */
  void
  reinit(const TriangulationManager<dim>    &triangulation_handler,
         const std::vector<FieldAttributes> &field_attributes);

  /**
   * @brief Getter function for the DoFHandlers (constant reference).
   */
  [[nodiscard]] const std::vector<const dealii::DoFHandler<dim> *> &
  get_dof_handlers(const std::set<unsigned int> &field_indices,
                   unsigned int                  relative_level = 0) const;

  /**
   * @brief Getter function for the DoFHandler (reference).
   */
  [[nodiscard]] const dealii::DoFHandler<dim> &
  get_dof_handler(Types::Index index, unsigned int relative_level = 0) const
  {
    return *dof_handlers[index][relative_level];
  }

  /**
   * @brief Get the total DoFs excluding multigrid DoFs.
   */
  [[nodiscard]] dealii::types::global_dof_index
  get_total_dofs() const
  {
    dealii::types::global_dof_index n_dofs = 0;
    for (const auto &dof_handler_set : dof_handlers)
      {
        n_dofs += dof_handler_set[0]->n_dofs();
      }
    return n_dofs;
  }

private:
  /**
   * @brief Collection of the triangulation DoFs. The number of DoFHandlers should be
   * equal to or less than the number of fields. Technically, there's a small
   * optimization we can use when multiple fields have the same constraints and
   * quadrature rule, allowing us to share the same DoFHandler. An example of this might
   * be grain growth.
   * Outer vector is indexed by field index. Inner vector is indexed by relative mg level.
   */
  std::vector<std::vector<std::shared_ptr<dealii::DoFHandler<dim>>>> dof_handlers;

#ifdef ADDITIONAL_OPTIMIZATIONS
  std::vector<unsigned int> smallest_index_with_identical_constraints;
#endif

  // TODO (fractalsbyx): move somewhere else. Maybe make an array?
  static const std::map<FieldInfo::TensorRank, dealii::FESystem<dim>> fe_systems {
    {FieldInfo::TensorRank::Scalar, dealii::FESystem<dim>(dealii::FE_Q<dim>(1), 1)  },
    {FieldInfo::TensorRank::Vector, dealii::FESystem<dim>(dealii::FE_Q<dim>(1), dim)}
  };
};

PRISMS_PF_END_NAMESPACE
