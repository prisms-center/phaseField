/*
 * nucleus.h
 *
 *  Created on: Mar 6, 2017
 *      Author: stephendewitt
 */

#ifndef APPLICATIONS__NUCLEATIONMODEL_NUCLEUS_H_
#define APPLICATIONS__NUCLEATIONMODEL_NUCLEUS_H_


// Structure representing each nucleus
struct nucleus{
    unsigned int index;
    dealii::Point<problemDIM> center;
    double radius;
    double seededTime, seedingTime;
    unsigned int seedingTimestep;
};



#endif /* APPLICATIONS__NUCLEATIONMODEL_NUCLEUS_H_ */
