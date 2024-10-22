#include "matrixFreePDE.h"

using namespace dealii;

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
  // Constructor
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

// Function to set the nucleation probability (in nucleation.h)
#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability([[maybe_unused]] variableValueContainer variable_value,
                           [[maybe_unused]] double                 dV) const override;
#endif

  // ================================================================
  // Methods specific to this subclass
  // ================================================================

  // Method that caps the value of the order parameter and the domain parameter
  void
  capFields(VectorizedArray<double> &ncp,
            VectorizedArray<double> &psicp,
            VectorizedArray<double>  n,
            VectorizedArray<double>  psi) const;

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  double WV       = userInputs.get_model_constant_double("WV");
  double gammaV   = userInputs.get_model_constant_double("gammaV");
  double icorrV   = userInputs.get_model_constant_double("icorrV");
  double VMV      = userInputs.get_model_constant_double("VMV");
  double zMV      = userInputs.get_model_constant_double("zMV");
  double zPV      = userInputs.get_model_constant_double("zPV");
  double znV      = userInputs.get_model_constant_double("znV");
  double DMV      = userInputs.get_model_constant_double("DMV");
  double DPV      = userInputs.get_model_constant_double("DPV");
  double DnV      = userInputs.get_model_constant_double("DnV");
  double cMsatV   = userInputs.get_model_constant_double("cMsatV");
  double epssqV   = userInputs.get_model_constant_double("epssqV");
  double EcorrV   = userInputs.get_model_constant_double("EcorrV");
  double VsV      = userInputs.get_model_constant_double("VsV");
  double betaV    = userInputs.get_model_constant_double("betaV");
  double TV       = userInputs.get_model_constant_double("TV");
  double rad0     = userInputs.get_model_constant_double("rad0");
  double lthresh  = userInputs.get_model_constant_double("lthresh");
  double imax_max = userInputs.get_model_constant_double("imax_max");
  double imax_min = userInputs.get_model_constant_double("imax_min");

  double FarC      = 96485.33289;
  double RV        = 8.314;
  double expconstV = zMV * (1.0 - betaV) * FarC / (RV * TV);
  double MconstV   = 2.0 * (VMV / (zMV * FarC)) * std::sqrt(2.0 * epssqV / WV);
  double deltaV    = std::sqrt(2.0 * epssqV / WV);

  // ================================================================
};
