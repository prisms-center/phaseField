#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
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

// Virtual method in MatrixFreePDE that we override if we need nucleation
#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability(variableValueContainer variable_value,
                           double                 dV,
                           dealii::Point<dim>     p,
                           unsigned int           variable_index) const;
#endif

  // ================================================================
  // Methods specific to this subclass
  // ================================================================

  // Method to place the nucleus and calculate the mobility modifier in
  // residualRHS
  void
  seedNucleus(const dealii::Point<dim, dealii::VectorizedArray<double>> &q_point_loc,
              dealii::VectorizedArray<double>                           &source_term,
              dealii::VectorizedArray<double>                           &gamma) const;

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
