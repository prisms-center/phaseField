#ifndef INCLUDE_ORDERPARAMETERREMAPPER_H_
#define INCLUDE_ORDERPARAMETERREMAPPER_H_

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/la_parallel_vector.h>

#include "SimplifiedGrainRepresentation.h"

#include <vector>

/**
 * This class uses information from the list of SimplifiedGrainRepresentation
 * objects to reassign grains across multiple solution fields.
 */
template <int dim>
class OrderParameterRemapper
{
public:
  /**
   * This method does the core work of the class to reassign grains across
   * solution vectors based on the list of SimplifiedGrainRepresentation
   * objects.
   */
  void
  remap(std::vector<SimplifiedGrainRepresentation<dim>> &grain_representations,
        std::vector<vectorType *>                       &solution_fields,
        dealii::DoFHandler<dim>                         &dof_handler,
        unsigned int                                     dofs_per_cell,
        double                                           buffer);

  /**
   * This method does the core work of the class to reassign grains across
   * solution vectors based on the list of SimplifiedGrainRepresentation
   * objects.
   */
  void
  remap_from_index_field(
    std::vector<SimplifiedGrainRepresentation<dim>> &grain_representations,
    const vectorType                                *grain_index_field,
    std::vector<vectorType *>                       &solution_fields,
    dealii::DoFHandler<dim>                         &dof_handler,
    unsigned int                                     dofs_per_cell,
    double                                           buffer);

protected:
};

#endif
