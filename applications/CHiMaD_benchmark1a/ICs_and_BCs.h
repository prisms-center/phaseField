// ===========================================================================
// FUNCTIONS FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim>
double InitialCondition<dim>::value (const Point<dim> &p, const unsigned int component) const
{
    double scalar_IC = 0;
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index.

    if (index == 0){
        double epsilon=0.01;
        double c0 = 0.5;
        double dx=userInputs.domain_size[0]/((double) userInputs.subdivisions[0])/std::pow(2.0,userInputs.refine_factor);
        scalar_IC = c0+epsilon*(std::cos(0.105*p[0])*std::cos(0.11*p[1])+
        std::pow(std::cos(0.13*p[0])*std::cos(0.087*p[1]),2.0)+
        std::cos(0.025*p[0]-0.15*p[1])*std::cos(0.07*p[0]-0.02*p[1]));

    }
    else {
        scalar_IC = 0.0;
    }

    // ---------------------------------------------------------------------
    return scalar_IC;
}

template <int dim>
void InitialConditionVec<dim>::vector_value (const Point<dim> &p,Vector<double> &vector_IC) const
{
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE FOR VECTOR FIELDS
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index.


    // ---------------------------------------------------------------------
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
