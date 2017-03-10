/*
 * nucleus.h
 *
 *  Created on: Mar 6, 2017
 *      Author: stephendewitt
 */

#ifndef APPLICATIONS__NUCLEATIONMODEL_NUCLEUS_H_
#define APPLICATIONS__NUCLEATIONMODEL_NUCLEUS_H_

#include "dealIIheaders.h"

// Structure representing each nucleus
template<int dim>
struct nucleus{
    unsigned int index;
    dealii::Point<dim> center;
    double radius;
    std::vector<double> semiaxes;
    double seededTime, seedingTime;
    unsigned int seedingTimestep;
};



#endif /* APPLICATIONS__NUCLEATIONMODEL_NUCLEUS_H_ */
