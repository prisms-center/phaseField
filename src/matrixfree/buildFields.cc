/*
 * buildFields.cc
 *
 *  Created on: Feb 22, 2017
 *      Author: stephendewitt
 */


// =====================================================================
// FUNCTION TO BUILD THE VECTOR OF FIELDS
// =====================================================================

template <int dim>
void MatrixFreePDE<dim>::buildFields(){
	// Build each of the fields in the system
	for (unsigned int i=0; i<num_var; i++){
		  if (var_type[i] == "SCALAR"){
			  if (var_eq_type[i] == "ELLIPTIC"){
				  this->fields.push_back(Field<problemDIM>(SCALAR, ELLIPTIC, var_name[i]));
			  }
			  else if (var_eq_type[i] == "PARABOLIC"){
				  this->fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, var_name[i]));
			  }
			  else{
				  // Need to change to throw an exception
				  std::cerr << "Error: Equation type must be ELLIPTIC or PARABOLIC " << std::endl;
			  }
		  }
		  else if (var_type[i] == "VECTOR"){
			  if (var_eq_type[i] == "ELLIPTIC"){
				  this->fields.push_back(Field<problemDIM>(VECTOR, ELLIPTIC, var_name[i]));
			  }
			  else if (var_eq_type[i] == "PARABOLIC"){
				  this->fields.push_back(Field<problemDIM>(VECTOR, PARABOLIC, var_name[i]));
			  }
			  else{
				  // Need to change to throw an exception
				  std::cerr << "Error: Variable type must be SCALAR or VECTOR " << std::endl;
			  }
		  }
	  }

}

