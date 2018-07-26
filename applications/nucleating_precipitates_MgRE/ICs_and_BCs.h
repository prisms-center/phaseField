// ===========================================================================
// FUNCTIONS FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim>
  double InitialCondition<dim>::value (const Point<dim> &p, const unsigned int component) const
  {
      double scalar_IC = 0;
      // --------------------------------------------------------------------------
      // ENTER THE INITIAL CONDITIONS HERE
      // --------------------------------------------------------------------------
      // Enter the function describing conditions for the fields at point "p".
      // Use "if" statements to set the initial condition for each variable
      // according to its variable index.

      // Initial condition parameters
      double c_matrix = 0.006;

      if (index==0){
          scalar_IC = c_matrix;
      }
      else if (index == 1){
          scalar_IC = 0.0;
      }
      else if (index==2||index==3){
          scalar_IC = 0.0;
      }
      

      // --------------------------------------------------------------------------
      return scalar_IC;
  }


//initial condition
template <int dim>
  void InitialConditionVec<dim>::vector_value (const Point<dim> &p,Vector<double> &vector_IC) const
  {
	  // --------------------------------------------------------------------------
	  // ENTER THE INITIAL CONDITIONS HERE FOR VECTOR FIELDS
	  // --------------------------------------------------------------------------
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  if (index==4){
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
