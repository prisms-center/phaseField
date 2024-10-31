#include "matrixFreePDE.h"

using namespace dealii;

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
  customPDE(userInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
    , userInputs(_userInputs) {};

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

private:
#include "typeDefs.h"

  const userInputParameters<dim> userInputs;

  // Function to set the RHS of the governing equations for explicit time
  // dependent equations (in equations.h)
  void
  explicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  // Function to set the RHS of the governing equations for all other equations
  // (in equations.h)
  void
  nonExplicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  // Function to set the LHS of the governing equations (in equations.h)
  void
  equationLHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

// Function to set postprocessing expressions (in postprocess.h)
#ifdef POSTPROCESS_FILE_EXISTS
  void
  postProcessedFields(
    [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
      &variable_list,
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &pp_variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;
#endif

// Virtual method in MatrixFreePDE that we override if we need nucleation
#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability([[maybe_unused]] variableValueContainer variable_value,
                           [[maybe_unused]] double                 dV,
                           [[maybe_unused]] Point<dim>             p,
                           [[maybe_unused]] unsigned int variable_index) const override;
#endif

  // ================================================================
  // Methods specific to this subclass
  // ================================================================

  // Method to place the nucleus and calculate the mobility modifier in
  // residualRHS
  void
  seedNucleus(const Point<dim, VectorizedArray<double>> &q_point_loc,
              VectorizedArray<double>                   &source_term,
              VectorizedArray<double>                   &gamma) const;

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  double c_avg     = userInputs.get_model_constant_double("c_avg");
  double McV       = userInputs.get_model_constant_double("McV");
  double MnV       = userInputs.get_model_constant_double("MnV");
  double KnV       = userInputs.get_model_constant_double("KnV");
  double W_barrier = userInputs.get_model_constant_double("W_barrier");
  double A0        = userInputs.get_model_constant_double("A0");
  double A2        = userInputs.get_model_constant_double("A2");
  double calmin    = userInputs.get_model_constant_double("calmin");
  double B0        = userInputs.get_model_constant_double("B0");
  double B2        = userInputs.get_model_constant_double("B2");
  double cbtmin    = userInputs.get_model_constant_double("cbtmin");

  double k1      = userInputs.get_model_constant_double("k1");
  double k2      = userInputs.get_model_constant_double("k2");
  double tau     = userInputs.get_model_constant_double("tau");
  double epsilon = userInputs.get_model_constant_double("epsilon");

  // Interface coefficient
  double interface_coeff = std::sqrt(2.0 * KnV / W_barrier);

  // ================================================================
};
