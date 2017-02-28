/*
 * varBCs.h
 *
 *  Created on: Feb 27, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_VARBCS_H_
#define INCLUDE_VARBCS_H_

// BC object declaration
template <int dim>
class varBCs
{
	public:
	//varBCs();
	std::vector<std::string> var_BC_type;
	std::vector<double> var_BC_val;
};



#endif /* INCLUDE_VARBCS_H_ */
