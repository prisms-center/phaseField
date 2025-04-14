// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>

#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
class userInputParameters;

template <int dim>
class triangulationHandler;

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
   * \brief Initialize the DoFHandlers
   */
  void
  init(const triangulationHandler<dim>                  &triangulation_handler,
       const std::map<fieldType, dealii::FESystem<dim>> &fe_system);

  /**
   * \brief Getter function for the DoFHandlers (constant reference).
   */
  [[nodiscard]] const std::vector<const dealii::DoFHandler<dim> *> &
  get_dof_handlers() const;

  /**
   * \brief Getter function for the multigrid DoFHandlers (constant reference).
   */
  [[nodiscard]] const std::map<unsigned int,
                               dealii::MGLevelObject<dealii::DoFHandler<dim>>> &
  get_mg_dof_handlers() const;

  /**
   * \brief Getter function for the DoFHandlers at a certain multigrid level
   * (constant reference).
   */
  [[nodiscard]] const std::vector<const dealii::DoFHandler<dim> *> &
  get_mg_dof_handlers(unsigned int level) const;

  /**
   * \brief Getter function for the DoFHandler at a certain field and multigrid level
   * (constant reference).
   */
  [[nodiscard]] const dealii::DoFHandler<dim> &
  get_mg_dof_handler(unsigned int index, unsigned int level) const;

private:
  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> *user_inputs;

  /**
   * \brief Collection of the triangulation DoFs. The number of DoFHandlers should be
   * equal to or less than the number of fields. Technically, there's a small
   * optimization we can use when multiple fields have the same constraints and
   * quadrature rule, allowing us to share the same DoFHandler. An example of this might
   * be grain growth.
   */
  std::map<unsigned int, std::unique_ptr<dealii::DoFHandler<dim>>> dof_handlers;

  /**
   * \brief Const copy of the dof_handlers.
   */
  std::vector<const dealii::DoFHandler<dim> *> const_dof_handlers;

  /**
   * \brief Whether we have multigrid.
   */
  bool has_multigrid = false;

  /**
   * \brief Collection of the triangulation DoFs for each multigrid level for all fields
   * that require it. Like before, we can share the same DoFHandler for multiple fields in
   * special cases.
   */
  std::map<unsigned int, dealii::MGLevelObject<dealii::DoFHandler<dim>>> mg_dof_handlers;

  /**
   * \brief Const copy of the mg_dof_handlers by multigrid level.
   */
  dealii::MGLevelObject<std::vector<const dealii::DoFHandler<dim> *>>
    const_mg_dof_handlers;
};

PRISMS_PF_END_NAMESPACE
