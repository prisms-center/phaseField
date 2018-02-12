// ===========================================================================
// FUNCTIONS FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim>
double InitialCondition<dim>::value (const dealii::Point<dim> &p, const unsigned int component) const
{
  double scalar_IC=0;
      // --------------------------------------------------------------------------
      // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
      // --------------------------------------------------------------------------
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  // Initial condition parameters
	  double x_denom = (2.0)*(2.0);
	  double y_denom = (3.0)*(3.0);
	  double z_denom = (12.0)*(12.0);
	  double initial_interface_coeff = 0.6;
	  double initial_radius = 1.0;
	  double c_matrix = 0.004;
      double c_precip = 0.25;

	  //set result equal to the structural order parameter initial condition
	  double r=0.0;
	  std::vector<double> ellipsoid_denoms;
	  ellipsoid_denoms.push_back(x_denom);
	  ellipsoid_denoms.push_back(y_denom);
	  ellipsoid_denoms.push_back(z_denom);
          
          dealii::Tensor<1,dim> center;
          center[0] = userInputs.domain_size[0]/2.0 + 1.234567;
          center[1] = userInputs.domain_size[1]/2.0 + 1.234567;

          if (dim == 3){
              center[2] = userInputs.domain_size[2]/2.0 + 1.234567;
          }
 /*
          std::vector<double> center_vec, semiaxes;
          center_vec.push_back(userInputs.domain_size[0]/2.0 + 1.234567);
          center_vec.push_back(userInputs.domain_size[1]/2.0 + 1.234567);
          if (dim == 3){
              center_vec.push_back(userInputs.domain_size[2]/2.0 + 1.234567);
          }
          dealii::Point<dim,double> center(center_vec);
   */
          std::vector<double> semiaxes;
          semiaxes.push_back(2.0);
          semiaxes.push_back(3.0);
          if (dim == 3){
              semiaxes.push_back(2.0);
          }

          std::vector<double> ellipsoid_rotation;
          ellipsoid_rotation.push_back(0.0);
          ellipsoid_rotation.push_back(0.0);
          ellipsoid_rotation.push_back(65.3);

          double degrees_to_rad = std::acos(0.0)/90.0;

          dealii::Tensor<2,dim,double> Rx, Ry, Rz;

          Rx[0][0] = 1.0;
          Rx[1][1] = std::cos(ellipsoid_rotation.at(0)*degrees_to_rad);

          Ry[0][0] = std::cos(ellipsoid_rotation.at(1)*degrees_to_rad);
          Ry[1][1] = 1.0;

          Rz[0][0] = std::cos(ellipsoid_rotation.at(2)*degrees_to_rad);
          Rz[1][0] = std::sin(ellipsoid_rotation.at(2)*degrees_to_rad);
          Rz[0][1] = -std::sin(ellipsoid_rotation.at(2)*degrees_to_rad);
          Rz[1][1] = std::cos(ellipsoid_rotation.at(2)*degrees_to_rad);

          if (dim == 3){
              Rx[1][2] = - std::sin(ellipsoid_rotation.at(0)*degrees_to_rad);
              Rx[2][1] = std::sin(ellipsoid_rotation.at(0)*degrees_to_rad);
              Rx[2][2] = std::cos(ellipsoid_rotation.at(0)*degrees_to_rad);

              Ry[0][2] = std::sin(ellipsoid_rotation.at(1)*degrees_to_rad);
              Ry[2][0] = -std::sin(ellipsoid_rotation.at(1)*degrees_to_rad);
              Ry[2][2] = std::cos(ellipsoid_rotation.at(1)*degrees_to_rad);

              Rz[2][2] = 1.0;
          }
          dealii::Tensor<2,dim,double> rotation_matrix;
          rotation_matrix = Rx*Ry*Rz; // Note: these are tensor multiplications

          // Calculate the distance function for a rotated ellipsoid
          double weighted_dist = 0.0;
          dealii::Tensor<1,dim,double> shortest_edist_tensor = center - p;
          for (unsigned int i=0; i<dim; i++){
              if (userInputs.BC_list[index].var_BC_type[2*i]==PERIODIC){
                  shortest_edist_tensor[i] = shortest_edist_tensor[i]-round(shortest_edist_tensor[i]/userInputs.domain_size[i])*userInputs.domain_size[i];
              }
          }
          shortest_edist_tensor = rotation_matrix * shortest_edist_tensor;
          for (unsigned int i=0; i<dim; i++){
             shortest_edist_tensor[i] /= semiaxes[i];
          }
          weighted_dist = shortest_edist_tensor.norm();          

          double avg_semiaxis = 0.0;
          for (unsigned int j=0; j<dim; j++){
              avg_semiaxis += semiaxes[j];
          }
          avg_semiaxis /= dim;

	  if (index==0){
	      scalar_IC = 0.5*(c_precip-c_matrix)*(1.0-std::tanh(avg_semiaxis*(weighted_dist-1.0)/initial_interface_coeff)) + c_matrix;	  
	  }
	  else if (index==1){
               scalar_IC = 0.5*(1.0-std::tanh(avg_semiaxis*(weighted_dist-1.0)/initial_interface_coeff));	
	  }

	  // --------------------------------------------------------------------------
	  return scalar_IC;
}

template <int dim>
void InitialConditionVec<dim>::vector_value (const dealii::Point<dim> &p, dealii::Vector<double> &vector_IC) const
{
	  // --------------------------------------------------------------------------
	  // ENTER THE INITIAL CONDITIONS HERE FOR VECTOR FIELDS
	  // --------------------------------------------------------------------------
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  if (index==2){
		  vector_IC(0) = 0.0;
		  vector_IC(1) = 0.0;
          if (dim == 3){
        	  vector_IC(2) = 0.0;
          }
	  }
	  // --------------------------------------------------------------------------
 }

 // ===========================================================================
 // FUNCTIONS FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
 // ===========================================================================

 template <int dim>
 double NonUniformDirichletBC<dim>::value (const dealii::Point<dim> &p, const unsigned int component) const
 {
     double scalar_BC=0;
     // --------------------------------------------------------------------------
     // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE FOR SCALAR FIELDS
     // --------------------------------------------------------------------------
     // Enter the function describing conditions for the fields at point "p".
     // Use "if" statements to set the boundary condition for each variable
     // according to its variable index. This function can be left blank if there
     // are no non-uniform Dirichlet boundary conditions.


     // -------------------------------------------------------------------------
     return scalar_BC;
 }

 template <int dim>
 void NonUniformDirichletBCVec<dim>::vector_value (const dealii::Point<dim> &p, dealii::Vector<double> &vector_BC) const
 {

     // --------------------------------------------------------------------------
     // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE FOR VECTOR FIELDS
     // --------------------------------------------------------------------------
     // Enter the function describing conditions for the fields at point "p".
     // Use "if" statements to set the boundary condition for each variable
     // according to its variable index. This function can be left blank if there
     // are no non-uniform Dirichlet boundary conditions.


     // -------------------------------------------------------------------------

 }
