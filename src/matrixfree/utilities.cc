//utility functions for the MatrixFreePDE class

#ifndef UTILITY_MATRIXFREE_H
#define UTILITY_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized


//return index of given field name if exists, else throw error
template <int dim>
unsigned int MatrixFreePDE<dim>::getFieldIndex(std::string _name) {
   for(typename std::vector<Field<dim> >::iterator it = fields.begin(); it != fields.end(); ++it){
     if (it->name.compare(_name)==0) return it->index;
   }
   pcout << "\nutilities.h: field '" << _name.c_str() << "' not initialized\n";
   exit(-1);
}

#endif
