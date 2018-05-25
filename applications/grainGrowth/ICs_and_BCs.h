// ===========================================================================
// FUNCTIONS FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim>
double InitialCondition<dim>::value (const dealii::Point<dim> &p, const unsigned int component) const
{
    double scalar_IC;
    // --------------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index.

    // The initial condition is a set of overlapping circles/spheres defined
    // by a hyperbolic tangent function. The center of each circle/sphere is
    // given by "center" and its radius is given by "radius".

    std::vector<dealii::Point<dim>> center;

    // The big grains
    {dealii::Point<dim> p(0.2, 0.15); center.push_back(p);}
    {dealii::Point<dim> p(0.25, 0.7); center.push_back(p);}
    {dealii::Point<dim> p(0.5, 0.5); center.push_back(p);}
    {dealii::Point<dim> p(0.6, 0.85); center.push_back(p);}
    {dealii::Point<dim> p(0.85, 0.35); center.push_back(p);}

    // The medium grains
    {dealii::Point<dim> p(0.08, 0.92); center.push_back(p);}
    {dealii::Point<dim> p(0.75, 0.6); center.push_back(p);}
    {dealii::Point<dim> p(0.75, 0.1); center.push_back(p);}
    {dealii::Point<dim> p(0.2, 0.45); center.push_back(p);}
    {dealii::Point<dim> p(0.85, 0.85); center.push_back(p);}

    // The small grains
    {dealii::Point<dim> p(0.55, 0.05); center.push_back(p);}
    {dealii::Point<dim> p(0.1, 0.35); center.push_back(p);}
    {dealii::Point<dim> p(0.95, 0.65); center.push_back(p);}
    {dealii::Point<dim> p(0.9, 0.15); center.push_back(p);}
    {dealii::Point<dim> p(0.45, 0.25); center.push_back(p);}


    //double center[15][3] = {{0.25,0.25,0},{0.25,0.5,0},{0.125,0.625,0},{0.375,0.625,0},{0.125,0.875,0},{0.375,0.875,0},{0.5,0.25,0},{0.5,0.5,0},{0.5,0.75,0},{0.75,0.25,0},{0.625,0.375,0},{0.875,0.375,0},{0.875,0.625,0},{0.875,0.625,0},{0.75,0.75,0}};

    std::vector<double> rad = {0.14, 0.14, 0.14, 0.14, 0.14, 0.08, 0.08, 0.08, 0.08, 0.08, 0.05, 0.05, 0.05, 0.05, 0.05};

    double dist = 0.0;
    scalar_IC = 0;

    for (unsigned int dir = 0; dir < dim; dir++){
        dist += (p[dir]-center[index][dir]*userInputs.domain_size[dir])*(p[dir]-center[index][dir]*userInputs.domain_size[dir]);
    }
    dist = std::sqrt(dist);


    scalar_IC += 0.5*(1.0-std::tanh((dist-rad[index]*userInputs.domain_size[0])/0.5));



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
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


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
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


    // -------------------------------------------------------------------------

}
