// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <prismspf/core/matrix_free_manager.h>
#include <prismspf/core/system_wide.h>

#include "prismspf/config.h"

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class holds information about a grain, including field index and the list
 * of vertices in that grain.
 */
template <unsigned int dim>
class GrainSet
{
public:
  /**
   * @brief Sets the grain index.
   */
  void
  set_grain_index(unsigned int _grain_index)
  {
    grain_index = _grain_index;
  };

  /**
   * @brief Gets the grain index.
   */
  [[nodiscard]] unsigned int
  get_grain_index() const
  {
    return grain_index;
  };

  /**
   * @brief Sets the order parameter index.
   */
  void
  set_order_parameter_index(unsigned int _op_index)
  {
    order_parameter_index = _op_index;
  };

  /**
   * @brief Gets the order parameter index.
   */
  [[nodiscard]] unsigned int
  get_order_parameter_index() const
  {
    return order_parameter_index;
  };

  /**
   * @brief Adds the vertices of a new element to the list.
   */
  void
  add_vertex_list(std::vector<dealii::Point<dim>> _vertices)
  {
    list_of_vertices.push_back(_vertices);
  };

  /**
   * @brief Gets the entire list of elements and the list of vertices per element.
   */
  std::vector<std::vector<dealii::Point<dim>>>
  get_vertex_list() const
  {
    return list_of_vertices;
  };

private:
  /**
   * @brief The grain index.
   */
  unsigned int grain_index;

  /**
   * @brief The variable index for the order parameter containing this grain.
   */
  unsigned int order_parameter_index;

  /**
   * @brief A vector of the elements in the grain containing a vector of the vertices
   * for each element.
   */
  std::vector<std::vector<dealii::Point<dim>>> list_of_vertices;
};

/**
 * This class uses a recursive flood filling algorithm to find connected bodies
 * in a field, given a threshold. The MPI communication methods are similar to
 * those in parallelNucleationList.
 */
template <unsigned int dim, unsigned int degree, typename number>
class FloodFiller
{
public:
  /**
   * @brief The primary external interface. This method takes in information about the
   * mesh/field and outputs a vector of GrainSet objects.
   */
  void
  calc_grain_sets(const dealii::DoFHandler<dim> &dof_handler,
                  const SolutionVector<number>  &solution_field,
                  double                         threshold_lower,
                  double                         threshold_upper,
                  int                            min_id,
                  unsigned int                   order_parameter_index,
                  std::vector<GrainSet<dim>>    &grain_sets);

private:
  /**
   * @brief The actual recursive flood fill method.
   */
  template <typename T>
  void
  recursive_flood_fill(T                             cell,
                       T                             cell_end,
                       const SolutionVector<number> &solution_field,
                       double                        threshold_lower,
                       double                        threshold_upper,
                       int                           min_id,
                       unsigned int                 &grain_index,
                       std::vector<GrainSet<dim>>   &grain_sets,
                       bool                         &grain_assigned);

  /**
   * @brief The method to merge the grain sets from all the processors.
   */
  void
  create_global_grain_set_list(std::vector<GrainSet<dim>> &grain_sets) const;

  /**
   * @brief Checks to see if grains found on different processors are parts of a larger
   * grain. If so, it merges the grain_sets entries.
   */
  void
  merge_split_grains(std::vector<GrainSet<dim>> &grain_sets) const;
};

PRISMS_PF_END_NAMESPACE