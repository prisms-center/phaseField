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
    double scalar_IC = 0;
    // =====================================================================
    // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
    // =====================================================================
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index.

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

    double mater_consts[1][2] = MaterialConstants;
    double E                  = mater_consts[0][0];
    double nu                 = mater_consts[0][1];
    double e0[dim][dim]       = {
      {0.01, 0,    0   },
      {0,    0.01, 0   },
      {0,    0,    0.01}
    };

    double a = 10.0;

    Point<dim> p_shifted;
    p_shifted[0] = p[0];
    p_shifted[1] = p[1];
    p_shifted[2] = p[2];

    double                 r = p_shifted.distance(Point<dim>(0, 0, 0));
    Tensor<2, dim, double> pT;
    Tensor<3, dim, double> S;

    for (unsigned int i = 0; i < dim; i++)
      {
        for (unsigned int j = 0; j < dim; j++)
          {
            pT[i][j] = E / (1.0 + nu) * (e0[i][j]);
            if (i == j)
              {
                pT[i][j] =
                  pT[i][j] + E / (1.0 + nu) *
                               (nu / (1.0 - 2.0 * nu) * (e0[0][0] + e0[1][1] + e0[2][2]));
              }
          }
      }

    if (r > a)
      {
        for (unsigned int i = 0; i < dim; i++)
          {
            double front_coeff = (1 + nu) * a * a * a / (2.0 * (1.0 - nu) * E);

            double term1 = 0;
            for (unsigned int k = 0; k < dim; k++)
              {
                term1 =
                  term1 + (3.0 * a * a - 5.0 * r * r) / (15.0 * r * r * r * r * r) *
                            (2.0 * pT[i][k] * p_shifted[k] + pT[k][k] * p_shifted[i]);
              }

            double term2 = 0;
            for (unsigned int j = 0; j < dim; j++)
              {
                for (unsigned int k = 0; k < dim; k++)
                  {
                    term2 = term2 + (r * r - a * a) / (r * r * r * r * r * r * r) *
                                      pT[j][k] * p_shifted[j] * p_shifted[i] *
                                      p_shifted[k];
                  }
              }

            double term3 = 0;
            for (unsigned int k = 0; k < dim; k++)
              {
                term3 =
                  term3 + 4.0 * (1 - nu) / (3.0 * r * r * r) * pT[i][k] * p_shifted[k];
              }

            vector_IC(i) = front_coeff * (term1 + term2 + term3);
          }
      }
    else
      {
        for (unsigned int i = 0; i < dim; i++)
          {
            for (unsigned int j = 0; j < dim; j++)
              {
                for (unsigned int k = 0; k < dim; k++)
                  {
                    for (unsigned int l = 0; l < dim; l++)
                      {
                        if ((i == j) & (i == k) & (i == l))
                          {
                            S[j][k][l] = (7.0 - 5.0 * nu) / (15.0 * (1.0 - nu));
                          }
                        else if ((i == j) & (k == l))
                          {
                            S[j][k][l] = (5.0 * nu - 1.0) / (15.0 * (1.0 - nu));
                          }
                        else if ((i == k) & (j == l))
                          {
                            S[j][k][l] = (4 - 5 * nu) / (15 * (1 - nu));
                          }
                        else if ((i == l) & (j == k))
                          {
                            S[j][k][l] = (4 - 5 * nu) / (15 * (1 - nu));
                          }
                        else
                          {
                            S[j][k][l] = 0.0;
                          }
                      }
                  }
              }

            vector_IC(i) = 0.0;

            for (unsigned int j = 0; j < dim; j++)
              {
                for (unsigned int k = 0; k < dim; k++)
                  {
                    for (unsigned int l = 0; l < dim; l++)
                      {
                        vector_IC(i) =
                          vector_IC(i) + S[j][k][l] * e0[k][l] * p_shifted[j];
                      }
                  }
              }
          }
      }

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

  inputBCs(4, 0, "ZERO_DERIVATIVE", 0.0);
  inputBCs(4, 1, "ZERO_DERIVATIVE", 0.0);
  inputBCs(4, 2, "ZERO_DERIVATIVE", 0.0);

  // =====================================================================
}
