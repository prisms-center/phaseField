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

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  double McV      = userInputs.get_model_constant_double("McV");
  double MnV      = userInputs.get_model_constant_double("MnV");
  double epsilonM = userInputs.get_model_constant_double("epsilonM");
  double delta2   = userInputs.get_model_constant_double("delta2");

  // ================================================================
};
