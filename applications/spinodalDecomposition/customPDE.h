#include "matrixFreePDE.h"

#include <random>

using namespace dealii;

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
  // Constructor
  customPDE(userInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
    , userInputs(_userInputs)
  {
    // Defining seed
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // Initializing distribution
    dist = distribution(0.0, 1.0);
    // Initializing random variable
    rng = engine(seed);
  };

  // Function to set the initial conditions (in ICs_and_BCs.h)
  void
  setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                      [[maybe_unused]] const unsigned int index,
                      [[maybe_unused]] double            &scalar_IC,
                      [[maybe_unused]] Vector<double>    &vector_IC) override;

  // Function to set the non-uniform Dirichlet boundary conditions (in
  // ICs_and_BCs.h)
  void
  setNonUniformDirichletBCs([[maybe_unused]] const Point<dim>  &p,
                            [[maybe_unused]] const unsigned int index,
                            [[maybe_unused]] const unsigned int direction,
                            [[maybe_unused]] const double       time,
                            [[maybe_unused]] double            &scalar_BC,
                            [[maybe_unused]] Vector<double>    &vector_BC) override;

  typedef std::mt19937_64                        engine;
  typedef std::uniform_real_distribution<double> distribution;

private:
#include "typeDefs.h"

  const userInputParameters<dim> userInputs;

  // Function to set the RHS of the governing equations for explicit time
  // dependent equations (in equations.h)
  void
  explicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                        &variable_list,
    [[maybe_unused]] Point<dim, VectorizedArray<double>> q_point_loc) const override;

  // Function to set the RHS of the governing equations for all other equations
  // (in equations.h)
  void
  nonExplicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                        &variable_list,
    [[maybe_unused]] Point<dim, VectorizedArray<double>> q_point_loc) const override;

  // Function to set the LHS of the governing equations (in equations.h)
  void
  equationLHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                        &variable_list,
    [[maybe_unused]] Point<dim, VectorizedArray<double>> q_point_loc) const override;

// Function to set postprocessing expressions (in postprocess.h)
#ifdef POSTPROCESS_FILE_EXISTS
  void
  postProcessedFields(
    [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
      &variable_list,
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &pp_variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc)
    const override;
#endif

// Function to set the nucleation probability (in nucleation.h)
#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability([[maybe_unused]] variableValueContainer variable_value,
                           [[maybe_unused]] double                 dV) const override;
#endif

  // ================================================================
  // Methods specific to this subclass
  // ================================================================

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  double McV         = userInputs.get_model_constant_double("McV");
  double KcV         = userInputs.get_model_constant_double("KcV");
  double WcV         = userInputs.get_model_constant_double("WcV");
  double c0          = userInputs.get_model_constant_double("c0");
  double icamplitude = userInputs.get_model_constant_double("icamplitude");

  // Declaring random number generator (Type std::mt19937_64)
  engine rng;
  // Declaring distribution (Type std::uniform_real_distribution<double>)
  distribution dist;

  // ================================================================
};
