#include "../../include/matrixFreePDE.h"

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

  // Function to set the non-uniform Dirichlet boundary conditions (in
  // ICs_and_BCs.h)
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

  // Function to set the RHS of the governing equations for explicit time
  // dependent equations (in equations.h)
  void
  explicitEquationRHS(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const;

  // Function to set the RHS of the governing equations for all other equations
  // (in equations.h)
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
#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability(variableValueContainer variable_value, double dV) const;
#endif

  // ================================================================
  // Methods specific to this subclass
  // ================================================================

  // Method that caps the value of the order parameter and the domain parameter
  void
  capFields(dealii::VectorizedArray<double> &ncp,
            dealii::VectorizedArray<double> &psicp,
            dealii::VectorizedArray<double>  n,
            dealii::VectorizedArray<double>  psi) const;

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
