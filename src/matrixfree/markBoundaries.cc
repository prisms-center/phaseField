//methods to mark boundaries 

#ifndef MARKBOUNDARIES_MATRIXFREE_H
#define MARKBOUNDARIES_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//methods to mark boundaries
//methods to mark boundaries
template <int dim>
void MatrixFreePDE<dim>::markBoundaries(){

	std::vector<double> domain_size;
	domain_size.push_back(spanX);
	domain_size.push_back(spanY);
	domain_size.push_back(spanZ);

	typename Triangulation<dim>::cell_iterator
	cell = MatrixFreePDE<dim>::triangulation.begin (),
	endc = MatrixFreePDE<dim>::triangulation.end();

	for (; cell!=endc; ++cell){

		// Mark all of the faces
		for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell;++face_number){
			for (unsigned int i=0; i<dim; i++){
				if ( std::fabs(cell->face(face_number)->center()(i) - (0)) < 1e-12 ){
					cell->face(face_number)->set_boundary_id (2*i);
				}
				else if (std::fabs(cell->face(face_number)->center()(i) - (domain_size[i])) < 1e-12){
					cell->face(face_number)->set_boundary_id (2*i+1);
				}

			}
		}
	}
}


#endif
