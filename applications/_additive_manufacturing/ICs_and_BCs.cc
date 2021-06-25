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

	  // The initial condition is a set of overlapping circles/spheres defined
	  // by a hyperbolic tangent function. The center of each circle/sphere is
	  // given by "center" and its radius is given by "radius".


          std::vector<dealii::Point<dim>> center;

          // The big grains
        //   {dealii::Point<dim> p(0.2, 0.15); center.push_back(p);}
        //   {dealii::Point<dim> p(0.25, 0.7); center.push_back(p);}
        //   {dealii::Point<dim> p(0.5, 0.5); center.push_back(p);}
        //   {dealii::Point<dim> p(0.6, 0.85); center.push_back(p);}
        //   {dealii::Point<dim> p(0.85, 0.35); center.push_back(p);}

        //   {dealii::Point<dim> p(0.08, 0.92); center.push_back(p);}
        //   {dealii::Point<dim> p(0.75, 0.6); center.push_back(p);}
        //   {dealii::Point<dim> p(0.75, 0.1); center.push_back(p);}
        //   {dealii::Point<dim> p(0.2, 0.45); center.push_back(p);}

          {dealii::Point<dim> p(0.0, 0.15); center.push_back(p);}
          {dealii::Point<dim> p(0.0, 0.7); center.push_back(p);}
          {dealii::Point<dim> p(0.0, 0.5); center.push_back(p);}
          {dealii::Point<dim> p(0.0, 0.85); center.push_back(p);}
          {dealii::Point<dim> p(0.0, 0.35); center.push_back(p);}

          {dealii::Point<dim> p(0.0, 0.92); center.push_back(p);}
          {dealii::Point<dim> p(0.0, 0.6); center.push_back(p);}
          {dealii::Point<dim> p(0.0, 0.1); center.push_back(p);}
          {dealii::Point<dim> p(0.0, 0.45); center.push_back(p);}

double init_factor = 1.0e6;


          std::vector<double> rad = {0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.025, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05};


    	  double dist = 0.0;
    	  scalar_IC = 0;
if (index==0){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[0][dir]*userInputs.domain_size[dir])*(p[dir]-center[0][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[0]*userInputs.domain_size[0])/0.5));
    }
dist=0.0;
    if (index==1){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[1][dir]*userInputs.domain_size[dir])*(p[dir]-center[1][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[1]*userInputs.domain_size[0])/0.5));
    }

    dist=0.0;
    if (index==2){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[2][dir]*userInputs.domain_size[dir])*(p[dir]-center[2][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[2]*userInputs.domain_size[0])/0.5));
    }

        dist=0.0;
    if (index==3){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[3][dir]*userInputs.domain_size[dir])*(p[dir]-center[3][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[3]*userInputs.domain_size[0])/0.5));
    }

        dist=0.0;
    if (index==4){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[4][dir]*userInputs.domain_size[dir])*(p[dir]-center[4][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[4]*userInputs.domain_size[0])/0.5));
    }
    dist=0.0;
    if (index==5){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[5][dir]*userInputs.domain_size[dir])*(p[dir]-center[5][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[5]*userInputs.domain_size[0])/0.5));
    }
        dist=0.0;
    if (index==0){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[6][dir]*userInputs.domain_size[dir])*(p[dir]-center[6][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[6]*userInputs.domain_size[0])/0.5));
    }

            dist=0.0;
    if (index==0){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[7][dir]*userInputs.domain_size[dir])*(p[dir]-center[7][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[7]*userInputs.domain_size[0])/0.5));
    }

    //             dist=0.0;
    // if (index==0){
    // 	  for (unsigned int dir = 0; dir < dim; dir++){
    // 		  dist += (p[dir]-center[8][dir]*userInputs.domain_size[dir])*(p[dir]-center[8][dir]*userInputs.domain_size[dir]);
    // 	  }
    // 	  dist = std::sqrt(dist);

    // 	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[8]*userInputs.domain_size[0])/0.5));
    // }









dist=0.0;
    if (index==6){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[0][dir]*userInputs.domain_size[dir])*(p[dir]-center[0][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[0]*userInputs.domain_size[0])/0.5));
    }

    dist=0.0;
    if (index==6){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[1][dir]*userInputs.domain_size[dir])*(p[dir]-center[1][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[1]*userInputs.domain_size[0])/0.5));
    }

        dist=0.0;
    if (index==6){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[2][dir]*userInputs.domain_size[dir])*(p[dir]-center[2][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[2]*userInputs.domain_size[0])/0.5));
    }

            dist=0.0;
    if (index==6){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[3][dir]*userInputs.domain_size[dir])*(p[dir]-center[3][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[3]*userInputs.domain_size[0])/0.5));
    }

            dist=0.0;
    if (index==6){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[4][dir]*userInputs.domain_size[dir])*(p[dir]-center[4][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[4]*userInputs.domain_size[0])/0.5));
    }

            dist=0.0;
    if (index==6){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[5][dir]*userInputs.domain_size[dir])*(p[dir]-center[5][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[5]*userInputs.domain_size[0])/0.5));
    }

                dist=0.0;
    if (index==6){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[6][dir]*userInputs.domain_size[dir])*(p[dir]-center[6][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[6]*userInputs.domain_size[0])/0.5));
    }

                    dist=0.0;
    if (index==6){
    	  for (unsigned int dir = 0; dir < dim; dir++){
    		  dist += (p[dir]-center[7][dir]*userInputs.domain_size[dir])*(p[dir]-center[7][dir]*userInputs.domain_size[dir]);
    	  }
    	  dist = std::sqrt(dist);

    	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[7]*userInputs.domain_size[0])/0.5));
    }

    //                 dist=0.0;
    // if (index==6){
    // 	  for (unsigned int dir = 0; dir < dim; dir++){
    // 		  dist += (p[dir]-center[8][dir]*userInputs.domain_size[dir])*(p[dir]-center[8][dir]*userInputs.domain_size[dir]);
    // 	  }
    // 	  dist = std::sqrt(dist);

    // 	  scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad[8]*userInputs.domain_size[0])/0.5));
    // }


	  // --------------------------------------------------------------------------



	  if (index==7){
		//   if(p[0]<userInputs.domain_size[0]*0.025){
		//   scalar_IC=1545.0;
		//   }
		//   else{

			scalar_IC=1551.0;

		//   }
		  }


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
