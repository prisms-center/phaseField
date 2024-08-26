/*
 * vectorLoad.h
 *
 *  Created on: Feb 22, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_VECTORLOAD_H_
#define INCLUDE_VECTORLOAD_H_

void
vectorLoad(bool in[], int array_size, std::vector<bool> &out);
void
vectorLoad(double in[], int array_size, std::vector<double> &out);
void
vectorLoad(int in[], int array_size, std::vector<int> &out);
void
vectorLoad(unsigned int in[], int array_size, std::vector<unsigned int> &out);
void
vectorLoad(std::string in[], int array_size, std::vector<std::string> &out);

#endif /* INCLUDE_VECTORLOAD_H_ */
