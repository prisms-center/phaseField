#ifndef INCLUDE_ORDERPARAMETERREMAPPER_H_
#define INCLUDE_ORDERPARAMETERREMAPPER_H_

#include "dealIIheaders.h"

/**
* This class uses information from the list of SimplifiedGrainRepresentation objects to reassign grains across multiple solution fields. (Currently this class is a one-function stub, it may be reassigned as part of a different class if it stays this way.)
*/
template <int dim>
class OrderParameterRemapper
{
public:

    /**
    * This method does the core work of the class to reassign grains across solution vectors based on the list of SimplifiedGrainRepresentation objects.
    */
    void remap(std::vector<SimplifiedGrainRepresentation<dim>> & grain_representations, std::vector<vectorType*> & solution_fields, dealii::DoFHandler<dim> &dof_handler, unsigned int dofs_per_cell);


protected:

};

#endif
