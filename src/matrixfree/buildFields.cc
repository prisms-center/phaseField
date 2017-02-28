/*
 * buildFields.cc
 *
 *  Created on: Feb 22, 2017
 *      Author: stephendewitt
 */


// =====================================================================
// FUNCTION TO BUILD THE VECTOR OF FIELDS
// =====================================================================

template <int dim, int degree>
void MatrixFreePDE<dim,degree>::buildFields(){
	// Build each of the fields in the system
	for (unsigned int i=0; i<userInputs.number_of_variables; i++){
		  if (userInputs.var_type[i] == "SCALAR"){
			  if (userInputs.var_eq_type[i] == "ELLIPTIC"){
				  this->fields.push_back(Field<dim>(SCALAR, ELLIPTIC, userInputs.var_name[i]));
			  }
			  else if (userInputs.var_eq_type[i] == "PARABOLIC"){
				  this->fields.push_back(Field<dim>(SCALAR, PARABOLIC, userInputs.var_name[i]));
			  }
			  else{
				  // Need to change to throw an exception
				  std::cerr << "Error: Equation type must be ELLIPTIC or PARABOLIC " << std::endl;
			  }
		  }
		  else if (userInputs.var_type[i] == "VECTOR"){
			  if (userInputs.var_eq_type[i] == "ELLIPTIC"){
				  this->fields.push_back(Field<dim>(VECTOR, ELLIPTIC, userInputs.var_name[i]));
			  }
			  else if (userInputs.var_eq_type[i] == "PARABOLIC"){
				  this->fields.push_back(Field<dim>(VECTOR, PARABOLIC, userInputs.var_name[i]));
			  }
			  else{
				  // Need to change to throw an exception
				  std::cerr << "Error: Variable type must be SCALAR or VECTOR " << std::endl;
			  }
		  }
	  }

}

