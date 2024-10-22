#include <deal.II/base/exceptions.h>

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

  void
  makeTriangulation(parallel::distributed::Triangulation<dim> &) const override;

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  scalarvalueType c_alpha = constV(0.3);
  scalarvalueType c_beta  = constV(0.7);
  scalarvalueType rho_s   = constV(5.0);
  scalarvalueType kappa   = constV(2.0);
  scalarvalueType M       = constV(5.0);

  // ================================================================
};

#include <deal.II/grid/grid_generator.h>

template <int dim, int degree>
void
customPDE<dim, degree>::makeTriangulation(
  parallel::distributed::Triangulation<dim> &tria) const
{
  parallel::distributed::Triangulation<dim> tria_horizontal_box(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> tria_vertical_box(MPI_COMM_WORLD);

  // Box dimensions
  double a = 100;
  double c = 20;

  // Check that dimensions match the benchmark
  AssertThrow(dim == 2, ExcMessage("CHiMaD Benchmark 1c should only be run in 2D."));

  // Overriding user specified subdivisions. This way we don't have to move points on the
  // mesh to match the two triangulations.
  AssertThrow(userInputs.subdivisions[0] == userInputs.subdivisions[1],
              ExcMessage("Subdivisions are automatically determined in this "
                         "application. Please make sure the x and y subdivisions match "
                         "in order to change the refinement."));
  std::vector<unsigned int> horizontal_subdivisions(dim, userInputs.subdivisions[0]);
  std::vector<unsigned int> vertical_subdivisions(dim, userInputs.subdivisions[0]);
  for (unsigned int i = 0; i < dim; i++)
    {
      if (i == 0)
        {
          horizontal_subdivisions[i] = 5 * userInputs.subdivisions[0];
        }
      else if (i == 1)
        {
          vertical_subdivisions[i] = 5 * userInputs.subdivisions[0];
        }
    }

  // Create bounding points for each part of the triangulation
  Point<dim> horizontal_origin;
  Point<dim> horizontal_corner;
  Point<dim> vertical_origin;
  Point<dim> vertical_corner;

  if (dim == 2)
    {
      horizontal_origin = Point<dim>(0.0, a);
      horizontal_corner = Point<dim>(a, a + c);
      vertical_origin   = Point<dim>(0.5 * (a - c), 0.0);
      vertical_corner   = Point<dim>(0.5 * (a + c), a);
    }

  GridGenerator::subdivided_hyper_rectangle(tria_horizontal_box,
                                            horizontal_subdivisions,
                                            horizontal_origin,
                                            horizontal_corner);
  GridGenerator::subdivided_hyper_rectangle(tria_vertical_box,
                                            vertical_subdivisions,
                                            vertical_origin,
                                            vertical_corner);

  // Merge the two triangulations
  GridGenerator::merge_triangulations(tria_horizontal_box, tria_vertical_box, tria);

  // Mark the boundaries
  for (const auto &cell : tria.active_cell_iterators())
    {
      // Mark all of the faces on the boundary with a boundary id of 0. This reduces the
      // complexity of the code at the cost of flexibility in boundary conditions. For the
      // benchmark case, we don't care about flexibility. If you plan to use this code to
      // create your own triangulation, modify this section accordingly.
      for (unsigned int face_number = 0; face_number < GeometryInfo<dim>::faces_per_cell;
           ++face_number)
        {
          const auto &face = cell->face(face_number);

          if (face->at_boundary())
            {
              face->set_boundary_id(0);
            }
        }
    }
}