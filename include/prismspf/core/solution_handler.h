// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

#include <map>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
class matrixfreeHandler;

struct variableAttributes;

/**
 * \brief Class that manages solution initialization and swapping with old solutions.
 */
template <unsigned int dim>
class solutionHandler
{
public:
  using VectorType   = dealii::LinearAlgebra::distributed::Vector<double>;
  using MGVectorType = dealii::LinearAlgebra::distributed::Vector<float>;

  /**
   * \brief Constructor.
   */
  solutionHandler(const std::map<unsigned int, variableAttributes> &_attributes_list,
                  const MGInfo<dim>                                &_mg_info);

  /**
   * \brief Get the solution vector set. This contains all the normal fields and is
   * typically used for output.
   *
   * TODO (landinjm): Make const ptr?
   */
  [[nodiscard]] std::map<unsigned int, VectorType *>
  get_solution_vector() const;

  /**
   * \brief Get a solution vector of a given field index and dependency type.
   *
   * TODO (landinjm): Make const ptr?
   */
  [[nodiscard]] VectorType *
  get_solution_vector(unsigned int index, DependencyType dependency_type) const;

  /**
   * \brief Get the "new" solution vector set.
   *
   * TODO (landinjm): Make const ptr?
   */
  [[nodiscard]] std::map<unsigned int, VectorType *>
  get_new_solution_vector() const;

  /**
   * \brief Get the "new" solution vector of a given field index.
   *
   * TODO (landinjm): Make const ptr?
   */
  [[nodiscard]] VectorType *
  get_new_solution_vector(unsigned int index) const;

  /**
   * \brief Get the mg solution vector set at a given level.
   */
  [[nodiscard]] std::vector<MGVectorType *>
  get_mg_solution_vector(unsigned int level) const;

  /**
   * \brief Get the mg solution vector set at a given level and index;
   */
  [[nodiscard]] MGVectorType *
  get_mg_solution_vector(unsigned int level, unsigned int index) const;

  /**
   * \brief Initialize the solution set.
   */
  void
  init(matrixfreeHandler<dim, double> &matrix_free_handler);

  /**
   * \brief Initialize the multigrid solution set.
   */
  void
  mg_init(
    const dealii::MGLevelObject<matrixfreeHandler<dim, float>> &mg_matrix_free_handler);

  /**
   * \brief Update the ghost values.
   *
   * TODO (landinjm): Fix so this isn't as wasteful in updating ghost values for all
   * solution vectors.
   */
  void
  update_ghosts() const;

  /**
   * \brief Apply the given constraints to a solution vector of a given field index.
   *
   * Note this applies constraints for all dependencyTypes of the given index.
   */
  void
  apply_constraints(unsigned int                             index,
                    const dealii::AffineConstraints<double> &constraints);

  /**
   * \brief Apply intial condition to the old fields. For now, this simply copies the
   * values in the normal field to the old.
   *
   * TODO (landinjm): What should we do for the initial condition of old fields.
   */
  void
  apply_initial_condition_for_old_fields();

  /**
   * \brief Update the `solution_set` with the `new_solution_set`. This has different
   * variants on which solutions to swap based on the FieldSolveType.
   */
  void
  update(const FieldSolveType &field_solve_type, const unsigned int &variable_index = 0);

private:
  /**
   * \brief The attribute list of the relevant variables.
   */
  const std::map<unsigned int, variableAttributes> *attributes_list;

  /**
   * \brief Whether multigrid has been enabled.
   */
  bool has_multigrid = false;

  /**
   * \brief Global minimum level for multigrid.
   */
  unsigned int global_min_level;

  /**
   * \brief Multigrid information.
   */
  const MGInfo<dim> *mg_info;

  /**
   * \brief The collection of solution vector at the current timestep. This includes
   * current values and old values.
   */
  std::map<std::pair<unsigned int, DependencyType>, std::unique_ptr<VectorType>>
    solution_set;

  /**
   * \brief The collection of new solution vectors at the current timestep. This is the
   * dst vector that is filled in the cell_loop. Unlike before, this only include the
   * current values which get updated in the `solution_set`.
   */
  std::map<unsigned int, std::unique_ptr<VectorType>> new_solution_set;

  /**
   * \brief The collection of solution vectors at the current timestep for the multigrid
   * hierarchy.
   */
  std::vector<std::vector<std::unique_ptr<MGVectorType>>> mg_solution_set;
};

PRISMS_PF_END_NAMESPACE
