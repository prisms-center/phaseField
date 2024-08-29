/*
 * varBCs.h
 *
 *  Created on: Feb 27, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_VARBCS_H_
#define INCLUDE_VARBCS_H_

// BC type enum
enum BC_type
{
  NATURAL,
  DIRICHLET,
  PERIODIC,
  NON_UNIFORM_DIRICHLET,
  NEUMANN
};

// BC object declaration
template <int dim>
class varBCs
{
public:
  std::vector<BC_type> var_BC_type;
  std::vector<double>  var_BC_val;
};

#endif /* INCLUDE_VARBCS_H_ */
