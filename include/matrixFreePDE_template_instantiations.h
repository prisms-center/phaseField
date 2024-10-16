/*
 * matrixFreePDE_template_instantiations.h
 *
 *  Created on: Feb 28, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_MATRIXFREEPDE_TEMPLATE_INSTANTIATIONS_H_
#define INCLUDE_MATRIXFREEPDE_TEMPLATE_INSTANTIATIONS_H_

#ifndef MATRIXFREEPDE_TEMPLATE_INSTANTIATION
#  define MATRIXFREEPDE_TEMPLATE_INSTANTIATION

template class MatrixFreePDE<2, 1>;
template class MatrixFreePDE<3, 1>;

template class MatrixFreePDE<2, 2>;
template class MatrixFreePDE<3, 2>;

template class MatrixFreePDE<3, 3>;
template class MatrixFreePDE<2, 3>;

template class MatrixFreePDE<3, 4>;
template class MatrixFreePDE<2, 4>;

template class MatrixFreePDE<3, 5>;
template class MatrixFreePDE<2, 5>;

template class MatrixFreePDE<3, 6>;
template class MatrixFreePDE<2, 6>;
#endif

#endif /* INCLUDE_MATRIXFREEPDE_TEMPLATE_INSTANTIATIONS_H_ */
