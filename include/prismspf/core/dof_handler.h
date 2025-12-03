// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

template <unsigned int dim>
class TriangulationHandler;

template <unsigned int dim>
class MGInfo;

/**
 * @brief Class that manages the deal.II DoFHandlers
 */
template <unsigned int dim>
class DofHandler
{
public:
  /**
   * @brief Constructor.
   */
  DofHandler(const UserInputParameters<dim> &_user_inputs, const MGInfo<dim> &mg_info);

  /**
   * @brief Initialize the DoFHandlers
   */
  void
  init(const TriangulationHandler<dim> &triangulation_handler,
       const std::map<FieldInfo::TensorRank, dealii::FESystem<dim>> &fe_system,
       const MGInfo<dim>                                            &mg_info);

  /**
   * @brief Reinitialize the DoFHandlers
   */
  void
  reinit(const TriangulationHandler<dim> &triangulation_handler,
         const std::map<FieldInfo::TensorRank, dealii::FESystem<dim>> &fe_system,
         const MGInfo<dim>                                            &mg_info);

  /**
   * @brief Getter function for the DoFHandlers (constant reference).
   */
  [[nodiscard]] const std::vector<const dealii::DoFHandler<dim> *> &
  get_dof_handlers() const;

  /**
   * @brief Getter function for the DoFHandlers at a given multigrid level (constant
   * reference).
   */
  [[nodiscard]] const std::vector<const dealii::DoFHandler<dim> *> &
  get_mg_dof_handlers(unsigned int level) const;

  /**
   * @brief Getter function for the DoFHandler (reference).
   */
  [[nodiscard]] const dealii::DoFHandler<dim> &
  get_dof_handler(Types::Index index) const
  {
    Assert(dof_handlers.contains(index),
           dealii::ExcMessage(
             "Invalid index when trying to access a specific DoFHandler"));
    Assert(dof_handlers.at(index) != nullptr, dealii::ExcNotInitialized());
    return *dof_handlers.at(index);
  };

  /**
   * @brief Get the total DoFs excluding multigrid DoFs.
   */
  [[nodiscard]] dealii::types::global_dof_index
  get_total_dofs() const
  {
    dealii::types::global_dof_index n_dofs = 0;
    for (const auto &[index, dof_handler] : dof_handlers)
      {
        n_dofs += dof_handler->n_dofs();
      }
    return n_dofs;
  };

private:
  /**
   * @brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * @brief Collection of the triangulation DoFs. The number of DoFHandlers should be
   * equal to or less than the number of fields. Technically, there's a small
   * optimization we can use when multiple fields have the same constraints and
   * quadrature rule, allowing us to share the same DoFHandler. An example of this might
   * be grain growth.
   */
  std::map<unsigned int, std::unique_ptr<dealii::DoFHandler<dim>>> dof_handlers;

  /**
   * @brief Const copy of the dof_handlers.
   */
  std::vector<const dealii::DoFHandler<dim> *> const_dof_handlers;

  /**
   * @brief Whether we have multigrid.
   */
  bool has_multigrid = false;

  /**
   * @brief Global minimum level for multigrid.
   */
  unsigned int global_min_level = 0;

#ifdef ADDITIONAL_OPTIMIZATIONS
  /**
   * @brief A set of the degenerate field indices that are outside of the multigrid
   * dependencies.
   *
   * Since the degenerate field indices are solely determined by the FieldInfo::TensorRank
   * and the constraints, no info about any multigrid schemes is contained. This leads to
   * the edge case where the degenerate field index is not part of the multigrid LHS
   * dependencies. In this case, we have to keep track of that and initialize the
   * degenerate field index as if it were part of the multigrid dependencies.
   */
  std::set<Types::Index> degenerate_field_indices_outside_mg;
#endif

  /**
   * @brief Collection of the triangulation DoFs for each multigrid level for all fields
   * that require it. Like before, we can share the same DoFHandler for multiple fields in
   * special cases.
   */
  std::map<unsigned int, dealii::MGLevelObject<std::unique_ptr<dealii::DoFHandler<dim>>>>
    mg_dof_handlers;

  /**
   * @brief Const copy of the mg_dof_handlers.
   */
  std::vector<std::vector<const dealii::DoFHandler<dim> *>> const_mg_dof_handlers;
};

PRISMS_PF_END_NAMESPACE
