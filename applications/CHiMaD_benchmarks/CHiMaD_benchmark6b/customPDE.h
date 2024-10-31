#include <deal.II/base/tensor.h>
#include <deal.II/grid/manifold.h>

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

  // Virtual method in MatrixFreePDE
  void
  makeTriangulation(parallel::distributed::Triangulation<dim> &) const override;

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  double McV     = userInputs.get_model_constant_double("McV");
  double KcV     = userInputs.get_model_constant_double("KcV");
  double rho     = userInputs.get_model_constant_double("rho");
  double c_alpha = userInputs.get_model_constant_double("c_alpha");
  double c_beta  = userInputs.get_model_constant_double("c_beta");
  double k       = userInputs.get_model_constant_double("k");
  double epsilon = userInputs.get_model_constant_double("epsilon");

  // ================================================================
};

#include <deal.II/grid/grid_generator.h>

template <int dim, int degree>
void
customPDE<dim, degree>::makeTriangulation(
  parallel::distributed::Triangulation<dim> &tria) const
{
  parallel::distributed::Triangulation<dim> tria_box(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> tria_semicircle(MPI_COMM_WORLD);

  // Check that dimensions match the benchmark
  AssertThrow(dim == 2, ExcMessage("CHiMaD Benchmark 6b should only be run in 2D."));

  // Create bounding points for each part of the triangulation
  Point<dim> box_origin;
  Point<dim> box_corner;
  Point<dim> semicircle_origin;

  if (dim == 2)
    {
      box_corner = Point<dim>(userInputs.domain_size[0], userInputs.domain_size[1]);
      semicircle_origin =
        Point<dim>(userInputs.domain_size[0], userInputs.domain_size[1] / 2.0);
    }

  GridGenerator::subdivided_hyper_rectangle(tria_box,
                                            userInputs.subdivisions,
                                            box_origin,
                                            box_corner);

  GridGenerator::half_hyper_ball(tria_semicircle,
                                 semicircle_origin,
                                 userInputs.domain_size[1] / 2.0);

  // Find the two non-corner vertices on the right side of the rectangular mesh
  Point<dim> pt1;
  Point<dim> pt2;
  for (auto &cell : tria_box.active_cell_iterators())
    {
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        {
          Point<dim> &v = cell->vertex(i);
          if ((std::abs(v(0) - userInputs.domain_size[0]) < 1e-10) &&
              (v(1) > userInputs.domain_size[1] / 2.0) &&
              (v(1) < userInputs.domain_size[1] - 1.0e-10))
            {
              pt1 = v;
            }
          if ((std::abs(v(0) - userInputs.domain_size[0]) < 1e-10) &&
              (v(1) < userInputs.domain_size[1] / 2.0) && (v(1) > 1.0e-10))
            {
              pt2 = v;
            }
        }
    }

  // Move the vertices at the center of the half hyper ball so that they will
  // align with non-corner vertices on the right side of the rectangular mesh
  for (auto &cell : tria_semicircle.active_cell_iterators())
    {
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        {
          Point<dim> &v = cell->vertex(i);
          if ((std::abs(v(0) - userInputs.domain_size[0]) < 1e-10) &&
              (v(1) > userInputs.domain_size[1] / 2.0) &&
              (v(1) < userInputs.domain_size[1] - 1.0e-10))
            {
              v(1) = pt1(1);
            }
          if ((std::abs(v(0) - userInputs.domain_size[0]) < 1e-10) &&
              (v(1) < userInputs.domain_size[1] / 2.0) && (v(1) > 1.0e-10))
            {
              v(1) = pt2(1);
            }
        }
    }

  // Merge the rectangle and the semicircle
  GridGenerator::merge_triangulations(tria_box, tria_semicircle, tria);

  // Attach flat manifold to the entire domain
  tria.reset_all_manifolds();
  tria.set_manifold(0, FlatManifold<dim>());
  tria.set_all_manifold_ids(0);

  // Attach spherical manifold
  tria.set_manifold(8, SphericalManifold<dim>(semicircle_origin));

  // Set the 3 outer cells of semicircle to the spherical manifold
  for (const auto &cell : tria.active_cell_iterators())
    {
      const Point<dim> cell_center          = cell->center();
      const double     distance_from_center = cell_center.distance(semicircle_origin);

      if (cell_center[0] > userInputs.domain_size[0] + 1.0e-10 &&
          distance_from_center > 0.1 * userInputs.domain_size[1])
        {
          cell->set_all_manifold_ids(8);
        }
    }

  // Transfinite interpolation
  TransfiniteInterpolationManifold<dim> transfinite_manifold;
  transfinite_manifold.initialize(tria);
  tria.set_manifold(0, transfinite_manifold);

  // Mark the boundaries
  for (const auto &cell : tria.active_cell_iterators())
    {
      // Mark all of the faces
      for (unsigned int face_number = 0; face_number < GeometryInfo<dim>::faces_per_cell;
           ++face_number)
        {
          if (cell->face(face_number)->at_boundary())
            {
              for (unsigned int i = 0; i < dim; i++)
                {
                  if (i == 0)
                    {
                      if (std::fabs(cell->face(face_number)->center()(i) - (0)) < 1e-12)
                        {
                          cell->face(face_number)->set_boundary_id(2 * i);
                        }
                      else if (std::fabs(cell->face(face_number)->center()(i) >
                                         (userInputs.domain_size[i])))
                        {
                          cell->face(face_number)->set_boundary_id(2 * i + 1);
                        }
                    }
                  else
                    {
                      if (std::fabs(cell->face(face_number)->center()(i) - (0)) < 1e-12)
                        {
                          cell->face(face_number)->set_boundary_id(2 * i);
                        }
                      else if (std::fabs(cell->face(face_number)->center()(i) -
                                         (userInputs.domain_size[i])) < 1e-12)
                        {
                          cell->face(face_number)->set_boundary_id(2 * i + 1);
                        }
                    }
                }
            }
        }
    }
}
