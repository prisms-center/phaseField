// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/triangulation_manager.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Class that manages the deal.II DoFHandlers
 */
template <unsigned int dim, unsigned int degree>
class DoFManager
{
public:
  /**
   * @brief Constructor.
   */
  DoFManager() = default;

  /**
   * @brief Disable copying.
   */
  DoFManager(const DoFManager &)  = delete;
  DoFManager(const DoFManager &&) = delete;
  DoFManager
  operator=(const DoFManager &) = delete;
  DoFManager &
  operator=(DoFManager &&other) = delete;

  /**
   * @brief Destructor.
   */
  ~DoFManager() = default;

  /**
   * @brief Reinitialize the DoFHandlers
   * @param triangulation_manager The triangulation manager.
   * @param init_mg Whether to initialize the multigrid DoFHandlers.
   * @pre init() must have been called with the correct number of levels.
   * @note May invalidate existing DoFHandler pointers if number of levels changes.
   */
  void
  reinit(const TriangulationManager<dim> &triangulation_manager, bool with_mg = true);

  /**
   * @brief Reinitialize the DoFHandlers
   * @pre reinit() must have been called.
   */
  void
  reinit_mapping(const std::vector<FieldAttributes> &field_attributes);

  /**
   * @brief Getter function for all the DoFHandlers.
   * @pre reinit_mapping() must have been called.
   */
  [[nodiscard]] const std::vector<std::vector<const dealii::DoFHandler<dim> *>> &
  get_field_dof_handlers_levels() const;

  /**
   * @brief Getter function for all the DoFHandlers on a level.
   * @pre reinit_mapping() must have been called.
   */
  [[nodiscard]] const std::vector<const dealii::DoFHandler<dim> *> &
  get_field_dof_handlers(unsigned int relative_level = 0) const;

  /**
   * @brief Getter function for the DoFHandler (reference).
   * @pre reinit_mapping() must have been called.
   */
  [[nodiscard]] const dealii::DoFHandler<dim> &
  get_field_dof_handler(Types::Index field_index, unsigned int relative_level = 0) const;

  /**
   * @brief Getter the DoFHandlers for a block of fields.
   * @pre reinit_mapping() must have been called.
   */
  [[nodiscard]] std::vector<const dealii::DoFHandler<dim> *>
  get_block_dof_handlers(const std::set<unsigned int> &field_indices,
                         unsigned int                  relative_level = 0) const;

  /**
   * @brief Getter function for the scalar and vector DoFHandlers.
   */
  [[nodiscard]] const std::vector<std::array<dealii::DoFHandler<dim>, 2>> &
  get_dof_handlers_levels() const;

  /**
   * @brief Getter function for the scalar and vector DoFHandlers on a level.
   */
  [[nodiscard]] const std::array<dealii::DoFHandler<dim>, 2> &
  get_dof_handlers(unsigned int relative_level = 0) const;

  /**
   * @brief Getter function for a specific scalar or vector DoFHandler.
   */
  [[nodiscard]] const dealii::DoFHandler<dim> &
  get_dof_handler(const unsigned int &rank, unsigned int relative_level = 0) const;

  /**
   * @brief Get the total DoFs excluding multigrid DoFs.
   * @pre reinit_mapping() must have been called.
   */
  [[nodiscard]] dealii::types::global_dof_index
  get_total_dofs() const;

private:
  /**
   * @brief Pointers to the dof handlers for each field on every mg level.
   * Outer vector is indexed by relative mg level. Inner vector is indexed by field index.
   */
  std::vector<std::vector<const dealii::DoFHandler<dim> *>> field_dof_handlers;

  /**
   * @brief A scalar and a vector dof handler for each level
   */
  std::vector<std::array<dealii::DoFHandler<dim>, 2>> level_dof_handlers;
};

PRISMS_PF_END_NAMESPACE
