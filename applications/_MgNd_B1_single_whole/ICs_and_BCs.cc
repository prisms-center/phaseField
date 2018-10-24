// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index



	  double c_matrix = 0.004; //1.0e-6

      if (index==0){
          scalar_IC = c_matrix;
      }
      else if (index==1){
          // Initial condition parameters
          double initial_interface_coeff = 0.6; //0.6;
          double initial_radius = 1.0;
          double c_precip = 0.25;

          //set result equal to the structural order parameter initial condition
          double r=0.0;
            dealii::Tensor<1,dim> center;
            center[0] = userInputs.domain_size[0]/2.0;
            center[1] = userInputs.domain_size[1]/2.0;

            if (dim == 3){
                center[2] = userInputs.domain_size[2]/2.0;
            }

            std::vector<double> semiaxes;
            semiaxes.push_back(2.0);
            semiaxes.push_back(6.0);
            if (dim == 3){
                semiaxes.push_back(2.0);
            }

            std::vector<double> ellipsoid_rotation;
            ellipsoid_rotation.push_back(0.0);
            ellipsoid_rotation.push_back(0.0);
            ellipsoid_rotation.push_back(-5.3);

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

            scalar_IC = 0.5*(1.0-std::tanh(avg_semiaxis*(weighted_dist-1.0)/initial_interface_coeff));

      }
      else if (index < 4){
          scalar_IC = 0.0;
      }
      else {
          for (unsigned int d=0; d<dim; d++){
              vector_IC(d) = 0.0;
          }
      }

	  // --------------------------------------------------------------------------
  }



  // ===========================================================================
  // FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
  // ===========================================================================

  template <int dim, int degree>
  void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
  {
      // --------------------------------------------------------------------------
      // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
      // --------------------------------------------------------------------------
      // Enter the function describing conditions for the fields at point "p".
      // Use "if" statements to set the boundary condition for each variable
      // according to its variable index. This function can be left blank if there
      // are no non-uniform Dirichlet boundary conditions. For BCs that change in
      // time, you can access the current time through the variable "time". The
      // boundary index can be accessed via the variable "direction", which starts
      // at zero and uses the same order as the BC specification in parameters.in
      // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


      // -------------------------------------------------------------------------

  }
