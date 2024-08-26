// initial condition
template <int dim>
class InitialCondition : public Function<dim>
{
public:
  unsigned int   index;
  Vector<double> values;

  InitialCondition(const unsigned int _index)
    : Function<dim>(1)
    , index(_index)
  {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) + 1);
  }

  double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    double scalar_IC;
    // =====================================================================
    // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
    // =====================================================================
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index.

    // set result equal to the structural order paramter initial condition
    double dx = spanX / ((double) subdivisionsX) / std::pow(2.0, refineFactor);
    double dy = spanY / ((double) subdivisionsY) / std::pow(2.0, refineFactor);
    double dz = spanZ / ((double) subdivisionsZ) / std::pow(2.0, refineFactor);
    double r  = 0.0;

    double initial_radius          = 40.0;
    double initial_interface_coeff = 2.0;

#if problemDIM == 1
    r = p.operator()(0);
    return 0.5 * (1.0 - std::tanh((r - spanX / 16.0) / (0.1 * dx)));
#elif problemDIM == 2

    r = sqrt((p.operator()(0)) * (p.operator()(0)) / x_denom +
             (p.operator()(1)) * (p.operator()(1)) / y_denom);
    return 0.5 * (1.0 - std::tanh((r - initial_radius) / (initial_interface_coeff)));
// return 0.5*(1.0-std::tanh((r-spanX/16.0)/(3*dx)));
#elif problemDIM == 3

    // r=p.distance(Point<dim>(spanX/2.0,spanY/2.0,spanZ/2.0));
    // return 0.5*(1.0-std::tanh((r-spanX/8.0)/(3*dx)));

    r = sqrt((p.operator()(0) - spanX / 2.0) * (p.operator()(0) - spanX / 2.0) +
             (p.operator()(1) - spanY / 2.0) * (p.operator()(1) - spanY / 2.0) * 4.0 +
             (p.operator()(2) - spanZ / 2.0) * (p.operator()(2) - spanZ / 2.0));
    return 0.5 * (1.0 - std::tanh((r - initial_radius) / (initial_interface_coeff)));

// planar interface
// r=sqrt((p.operator()(2)-spanZ/2.0)*(p.operator()(2)-spanZ/2.0));
// return 0.5*(1.0-std::tanh((r)/(initial_interface_coeff)));
// return 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));
#endif
    return 0.0;

    // =====================================================================
    return scalar_IC;
  }
};

// initial condition
template <int dim>
class InitialConditionVec : public Function<dim>
{
public:
  unsigned int index;

  // Vector<double> values;
  InitialConditionVec(const unsigned int _index)
    : Function<dim>(dim)
    , index(_index)
  {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) + 1);
  }

  void
  vector_value(const Point<dim> &p, Vector<double> &vector_IC) const
  {
    // =====================================================================
    // ENTER THE INITIAL CONDITIONS HERE FOR VECTOR FIELDS
    // =====================================================================
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index.

    // =====================================================================
  }
};

// Sets the BCs for the problem variables
// "inputBCs" should be called for each component of each variable and should be
// in numerical order Four input arguments set the same BC on the entire
// boundary Two plus two times the number of dimensions inputs sets separate BCs
// on each face of the domain Inputs to "inputBCs": First input: variable number
// Second input: component number
// Third input: BC type (options are "ZERO_DERIVATIVE" and "DIRICHLET")
// Fourth input: BC value (ignored unless the BC type is "DIRICHLET")
// Odd inputs after the third: BC type
// Even inputs after the third: BC value
// Face numbering: starts at zero with the minimum of the first direction, one
// for the maximum of the first direction
//						two for the minimum of the second direction, etc.
template <int dim>
void
generalizedProblem<dim>::setBCs()
{
  // =====================================================================
  // ENTER THE BOUNDARY CONDITIONS HERE
  // =====================================================================

  inputBCs(0, 0, "ZERO_DERIVATIVE", 0);

  // =====================================================================
}
