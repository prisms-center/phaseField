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
              dealii::AlignedVector<dealii::VectorizedArray<double>>    &source_terms,
              dealii::VectorizedArray<double>                           &gamma) const;

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  double                 McV  = userInputs.get_model_constant_double("McV");
  double                 Mn1V = userInputs.get_model_constant_double("Mn1V");
  double                 Mn2V = userInputs.get_model_constant_double("Mn2V");
  double                 Mn3V = userInputs.get_model_constant_double("Mn3V");
  dealii::Tensor<2, dim> Kn1  = userInputs.get_model_constant_rank_2_tensor("Kn1");
  dealii::Tensor<2, dim> Kn2  = userInputs.get_model_constant_rank_2_tensor("Kn2");
  dealii::Tensor<2, dim> Kn3  = userInputs.get_model_constant_rank_2_tensor("Kn3");
  double                 W    = userInputs.get_model_constant_double("W");
  bool                   n_dependent_stiffness =
    userInputs.get_model_constant_bool("n_dependent_stiffness");
  dealii::Tensor<2, dim> sfts_const1 =
    userInputs.get_model_constant_rank_2_tensor("sfts_const1");
  dealii::Tensor<2, dim> sfts_const2 =
    userInputs.get_model_constant_rank_2_tensor("sfts_const2");
  dealii::Tensor<2, dim> sfts_const3 =
    userInputs.get_model_constant_rank_2_tensor("sfts_const3");
  double A2 = userInputs.get_model_constant_double("A2");
  double A1 = userInputs.get_model_constant_double("A1");
  double A0 = userInputs.get_model_constant_double("A0");
  double B2 = userInputs.get_model_constant_double("B2");
  double B1 = userInputs.get_model_constant_double("B1");
  double B0 = userInputs.get_model_constant_double("B0");

  const static unsigned int          CIJ_tensor_size = 2 * dim - 1 + dim / 3;
  dealii::Tensor<2, CIJ_tensor_size> CIJ_Mg =
    userInputs.get_model_constant_elasticity_tensor("CIJ_Mg");
  dealii::Tensor<2, CIJ_tensor_size> CIJ_Beta1 =
    userInputs.get_model_constant_elasticity_tensor("CIJ_Beta1");
  dealii::Tensor<2, CIJ_tensor_size> CIJ_Beta2 =
    userInputs.get_model_constant_elasticity_tensor("CIJ_Beta2");
  dealii::Tensor<2, CIJ_tensor_size> CIJ_Beta3 =
    userInputs.get_model_constant_elasticity_tensor("CIJ_Beta3");

  double interface_coeff = userInputs.get_model_constant_double("interface_coeff");
  double k1              = userInputs.get_model_constant_double("k1");
  double k2              = userInputs.get_model_constant_double("k2");
  double tau             = userInputs.get_model_constant_double("tau");

  double pi = 2.0 * std::acos(0.0);

  // ================================================================
};
