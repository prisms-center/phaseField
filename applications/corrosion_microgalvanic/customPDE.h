#include "../../include/matrixFreePDE.h"
#include <cmath>
#include <iostream>
using namespace std;

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
  setInitialCondition(const dealii::Point<dim> &p,
                      const unsigned int        index,
                      double                   &scalar_IC,
                      dealii::Vector<double>   &vector_IC);

  // Function to set the non-uniform Dirichlet
  // boundary conditions (in ICs_and_BCs.h)
  void
  setNonUniformDirichletBCs(const dealii::Point<dim> &p,
                            const unsigned int        index,
                            const unsigned int        direction,
                            const double              time,
                            double                   &scalar_BC,
                            dealii::Vector<double>   &vector_BC);

private:
#include "../../include/typeDefs.h"

  const userInputParameters<dim> userInputs;

  // Function to set the RHS of the governing equations
  // for explicit time dependent equations (in equations.h)

  void
  explicitEquationRHS(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const;

  // Function to set the RHS of the governing equations
  // for all other equations (in equations.h)
  void
  nonExplicitEquationRHS(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const;

  // Function to set the LHS of the governing equations (in equations.h)
  void
  equationLHS(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const;

  // Function to set postprocessing expressions (in postprocess.h)

#ifdef POSTPROCESS_FILE_EXISTS
  void
  postProcessedFields(
    const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &pp_variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<double>>        q_point_loc) const;
#endif

// Function to set the nucleation probability (in nucleation.h)
// Not useful for this application
#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability(variableValueContainer variable_value, double dV) const;
#endif

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  double VMV                = userInputs.get_model_constant_double("VMV");
  double zMV                = userInputs.get_model_constant_double("zMV");
  double epssqV             = userInputs.get_model_constant_double("epssqV");
  double EcorrAnodic        = userInputs.get_model_constant_double("EcorrAnodic");
  double EcorrCathodic      = userInputs.get_model_constant_double("EcorrCathodic");
  double AAnodic            = userInputs.get_model_constant_double("AAnodic");
  double ACathodic          = userInputs.get_model_constant_double("ACathodic");
  double VsV                = userInputs.get_model_constant_double("VsV");
  double lthresh            = userInputs.get_model_constant_double("lthresh");
  double gamma              = userInputs.get_model_constant_double("gamma");
  double kappa              = userInputs.get_model_constant_double("kappa");
  double i0Anodic           = userInputs.get_model_constant_double("i0Anodic");
  double i0Cathodic         = userInputs.get_model_constant_double("i0Cathodic");
  double iMax               = userInputs.get_model_constant_double("iMax");
  double cathodeThickness   = userInputs.get_model_constant_double("cathodeThickness");
  double tStepStartForV     = userInputs.get_model_constant_double("tStepStartForV");
  double initialGuessForPhi = userInputs.get_model_constant_double("guessValPhi");
  double FarC               = 96485.33289;        // Faraday constant [C/mol]
  double RV                 = 8.314;              // Gas constant [J/(mol K)]
  double deltaV             = sqrt(2.0 * epssqV); // Equilibrium interface half-width
  double MconstV = 2.0 * (VMV / (zMV * FarC)) * deltaV; // Constant prefactor for mobility

  // ================================================================
};
