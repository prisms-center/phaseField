#include "matrixFreePDE.h"

using namespace dealii;

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
  customPDE(userInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
    , userInputs(_userInputs)
  {
    c_dependent_misfit = false;
    for (unsigned int i = 0; i < dim; i++)
      {
        for (unsigned int j = 0; j < dim; j++)
          {
            if (std::abs(sfts_linear1[i][j]) > 1.0e-12)
              {
                c_dependent_misfit = true;
              }
          }
      }
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

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  double         McV  = userInputs.get_model_constant_double("McV");
  double         Mn1V = userInputs.get_model_constant_double("Mn1V");
  Tensor<2, dim> Kn1  = userInputs.get_model_constant_rank_2_tensor("Kn1");
  double         W    = userInputs.get_model_constant_double("W");
  bool           n_dependent_stiffness =
    userInputs.get_model_constant_bool("n_dependent_stiffness");
  Tensor<2, dim> sfts_linear1 =
    userInputs.get_model_constant_rank_2_tensor("sfts_linear1");
  Tensor<2, dim> sfts_const1 = userInputs.get_model_constant_rank_2_tensor("sfts_const1");

  double A2 = userInputs.get_model_constant_double("A2");
  double A1 = userInputs.get_model_constant_double("A1");
  double A0 = userInputs.get_model_constant_double("A0");
  double B2 = userInputs.get_model_constant_double("B2");
  double B1 = userInputs.get_model_constant_double("B1");
  double B0 = userInputs.get_model_constant_double("B0");

  const static unsigned int  CIJ_tensor_size = 2 * dim - 1 + dim / 3;
  Tensor<2, CIJ_tensor_size> CIJ_Mg =
    userInputs.get_model_constant_elasticity_tensor("CIJ_Mg");
  Tensor<2, CIJ_tensor_size> CIJ_Beta =
    userInputs.get_model_constant_elasticity_tensor("CIJ_Beta");

  bool c_dependent_misfit;

  // ================================================================
};
