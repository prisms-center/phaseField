#ifndef INCLUDE_ORDERPARAMETERREMAPPER_H_
#define INCLUDE_ORDERPARAMETERREMAPPER_H_

#include "dealIIheaders.h"

/**
* This class
*/
template <int dim>
class OrderParameterRemapper
{
public:
    void remap(std::vector<SimplifiedGrainRepresentation<dim>> & grain_representations, std::vector<vectorType*> & solution_fields, dealii::DoFHandler<dim> &dof_handler);


protected:

};

#endif
