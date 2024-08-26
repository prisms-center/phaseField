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

  // Method to refine the mesh
  void
  adaptiveRefineCriterion();

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

  double k1     = userInputs.get_model_constant_double("k1");
  double k2_b   = userInputs.get_model_constant_double("k2_b");
  double k2_gb  = userInputs.get_model_constant_double("k2_gb");
  double tau_b  = userInputs.get_model_constant_double("tau_b");
  double tau_gb = userInputs.get_model_constant_double("tau_gb");
  double wgb    = userInputs.get_model_constant_double("wgb");
  double gbll   = userInputs.get_model_constant_double("gbll");
  double gbrl   = userInputs.get_model_constant_double("gbrl");

  // Interface coefficient
  double interface_coeff = std::sqrt(2.0 * KnV / W_barrier);

  // ================================================================
};

// Special implementation of adaptive mesh criterion to make sure grain boundary
// region is adapted to the highest level
template <int dim, int degree>
void
customPDE<dim, degree>::adaptiveRefineCriterion()
{
  // Custom defined estimation criterion

  std::vector<std::vector<double>> errorOutV;

  QGaussLobatto<dim> quadrature(degree + 1);
  FEValues<dim>      fe_values(
    *(this->FESet[userInputs.refinement_criteria[0].variable_index]),
    quadrature,
    update_values | update_quadrature_points);
  const unsigned int              num_quad_points = quadrature.size();
  std::vector<dealii::Point<dim>> q_point_list(num_quad_points);

  std::vector<double> errorOut(num_quad_points);

  typename DoFHandler<dim>::active_cell_iterator
    cell = this->dofHandlersSet_nonconst[userInputs.refinement_criteria[0].variable_index]
             ->begin_active(),
    endc = this->dofHandlersSet_nonconst[userInputs.refinement_criteria[0].variable_index]
             ->end();

  typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell =
    this->triangulation.begin_active();

  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          for (unsigned int field_index = 0;
               field_index < userInputs.refinement_criteria.size();
               field_index++)
            {
              fe_values.get_function_values(
                *(this->solutionSet[userInputs.refinement_criteria[field_index]
                                      .variable_index]),
                errorOut);
              errorOutV.push_back(errorOut);
            }

          q_point_list = fe_values.get_quadrature_points();

          bool mark_refine = false;

          for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
            {
              for (unsigned int field_index = 0;
                   field_index < userInputs.refinement_criteria.size();
                   field_index++)
                {
                  bool cond_1 =
                    ((errorOutV[field_index][q_point] >
                      userInputs.refinement_criteria[field_index].value_lower_bound) &&
                     (errorOutV[field_index][q_point] <
                      userInputs.refinement_criteria[field_index].value_upper_bound));
                  bool cond_2 = (q_point_list[q_point](0) > gbll) &&
                                (q_point_list[q_point](0) < gbrl);
                  if (cond_1 || cond_2)
                    {
                      mark_refine = true;
                      break;
                    }
                }
            }

          errorOutV.clear();

          // limit the maximal and minimal refinement depth of the mesh
          unsigned int current_level = t_cell->level();

          if ((mark_refine && current_level < userInputs.max_refinement_level))
            {
              cell->set_refine_flag();
            }
          else if (!mark_refine && current_level > userInputs.min_refinement_level)
            {
              cell->set_coarsen_flag();
            }
        }
      ++t_cell;
    }
}
