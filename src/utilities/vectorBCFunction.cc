/*
 * vectorBCFunction.cc
 *
 *  Created on: Feb 22, 2017
 *      Author: stephendewitt
 */

template <int dim>
vectorBCFunction<dim>::vectorBCFunction(std::vector<double> input_values) : Function<dim>(dim), BC_values (input_values) {}

template <int dim>
void vectorBCFunction<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const {

	for (unsigned int i=0; i<dim; i++) {
		values(i) = BC_values[i];
	}
}

template <int dim>
void vectorBCFunction<dim>::vector_value_list (const std::vector<Point<dim> > &points, std::vector<Vector<double> > &value_list) const{
	const unsigned int n_points = points.size();
	for (unsigned int p=0; p<n_points; ++p)
		vectorBCFunction<dim>::vector_value(points[p],value_list[p]);
}


