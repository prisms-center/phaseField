// base class for matrix Free implementation of PDE's
#ifndef MATRIXFREEPDE_H
#define MATRIXFREEPDE_H

// dealii headers
#include <deal.II/base/quadrature.h>
#include <deal.II/base/timer.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR > 3)
#  include <deal.II/fe/mapping_fe.h>
#endif
#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/numerics/vector_tools.h>

// PRISMS headers
#include <core/fields.h>
#include <core/refinement/AdaptiveRefinement.h>
#include <core/userInputParameters.h>
#include <core/variableContainer.h>
#include <core/variableValueContainer.h>
#include <grains/SimplifiedGrainRepresentation.h>
#include <nucleation/nucleus.h>
#include <utilities/computeStress.h>

using namespace dealii;

// define data types
#ifndef scalarType
using scalarType = VectorizedArray<double>;
#endif
#ifndef vectorType
using vectorType = LinearAlgebra::distributed::Vector<double>;
#endif

// macro for constants
#define constV(a) make_vectorized_array(a)

/**
 * \brief This is the abstract base class for the matrix free implementation of parabolic
 * and elliptic BVP's, and supports MPI, threads and vectorization (Hybrid Parallel). This
 * class contains the parallel data structures, mesh (referred to as triangulation),
 * parallel degrees of freedom distribution, constraints, and general utility methods.
 *
 * All the physical models in this package inherit this base class.
 *
 * \tparam dim The number of dimensions in the problem.
 * \tparam degree The polynomial degree of the shape functions.
 */
template <int dim, int degree>
class MatrixFreePDE : public Subscriptor
{
public:
  /**
   * \brief Class contructor.
   */
  MatrixFreePDE(userInputParameters<dim>);

  /**
   * \brief Class destructor.
   */
  ~MatrixFreePDE() override;

  /**
   * \brief Initializes the mesh, degrees of freedom, constraints and data structures
   * using the user provided inputs in the application parameters file.
   */
  virtual void
  init();

  /**
   * \brief Create the mesh with the user provided domain sizes.
   *
   * \param tria Triangulation object. It must be empty prior to calling this function.
   */
  virtual void
  create_triangulation(parallel::distributed::Triangulation<dim> &tria) const;

  /**
   * \brief Initializes the data structures for enabling unit tests.
   *
   * \details This method initializes the MatrixFreePDE object with a fixed geometry,
   * discretization and other custom selected options specifically to help with unit
   * tests, and should not be called in any of the physical models.
   *
   * \param _fields Vector of PDE descriptions (e.g., scalar/vector) for each field.
   */
  void
  initForTests(std::vector<Field<dim>> _fields);

  /**
   * \brief This method implements the time stepping algorithm and invokes the
   * solveIncrement() method.
   */
  void
  solve();

  /**
   * \brief This method essentially converts the MatrixFreePDE object into a matrix object
   * which can be used with matrix free iterative solvers. Provides the A*x functionality
   * for solving the system of equations Ax=b.
   *
   * \param dst The destination vector.
   * \param src The source vector.
   */
  void
  vmult(vectorType &dst, const vectorType &src) const;

  /**
   * \brief Vector of all the physical fields in the problem. Fields are identified by
   * dimentionality (SCALAR/VECTOR), the kind of PDE (ELLIPTIC/PARABOLIC) used to compute
   * them and a character identifier  (e.g.: "c" for composition) which is used to write
   * the fields to the output files.
   */
  std::vector<Field<dim>> fields;

  /**
   * \brief Create the vector of all physical fields in the problem.
   */
  void
  buildFields();

  /**
   * \brief Parallel message stream.
   */
  ConditionalOStream pcout;

  /**
   * \brief Set the initial condition for all fields. This function is overriden in each
   * application.
   *
   * \param p The point at which the initial condition is evaluated.
   * \param index The index of the field being evaluated.
   * \param scalar_IC Return variable for scalar fields.
   * \param vector_IC Return variable for vector fields.
   */
  virtual void
  setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                      [[maybe_unused]] const unsigned int index,
                      [[maybe_unused]] double            &scalar_IC,
                      [[maybe_unused]] Vector<double>    &vector_IC) = 0;

  /**
   * \brief Set the spatially or temporally non-uniform boundary conditions. This function
   * is overriden in each application.
   *
   * \param p The point at which the boundary condition is evaluated.
   * \param index The index of the field being evaluated.
   * \param direction The face of the boundary condition.
   * \param time The time at which the boundary condition is evaluated.
   * \param scalar_BC Return variable for scalar fields.
   * \param vector_BC Return variable for vector fields.
   */
  virtual void
  setNonUniformDirichletBCs([[maybe_unused]] const Point<dim>  &p,
                            [[maybe_unused]] const unsigned int index,
                            [[maybe_unused]] const unsigned int direction,
                            [[maybe_unused]] const double       time,
                            [[maybe_unused]] double            &scalar_BC,
                            [[maybe_unused]] Vector<double>    &vector_BC) = 0;

protected:
  userInputParameters<dim> userInputs;

  unsigned int totalDOFs;

  // The attributes of the primary field variables and
  // the postprocessing field variables
  const AttributesList &var_attributes;
  const AttributesList &pp_attributes;

  // Elasticity matrix variables
  const static unsigned int CIJ_tensor_size = 2 * dim - 1 + dim / 3;

  // Method to reinitialize the mesh, degrees of freedom, constraints and data
  // structures when the mesh is adapted
  void
  reinit();

  /**
   * Method to reassign grains when multiple grains are stored in a single order
   * parameter.
   */
  void
  reassignGrains();

  std::vector<SimplifiedGrainRepresentation<dim>> simplified_grain_representations;

  /**
   * Method to solve each time increment of a time-dependent problem. For
   * time-independent problems this method is called only once. This method
   * solves for all the fields in a staggered manner (one after another) and
   * also invokes the corresponding solvers: Explicit solver for Parabolic
   * problems, Implicit (matrix-free) solver for Elliptic problems.
   */
  virtual void
  solveIncrement(bool skip_time_dependent);
  /* Method to write solution fields to vtu and pvtu (parallel) files.
   *
   * This method can be enabled/disabled by setting the flag writeOutput to
   * true/false. Also, the user can select how often the solution files are
   * written by setting the flag skipOutputSteps in the parameters file.
   */
  void
  outputResults();

  /*Parallel mesh object which holds information about the FE nodes, elements
   * and parallel domain decomposition
   */
  parallel::distributed::Triangulation<dim> triangulation;
  /*A vector of finite element objects used in a model. For problems with only
   *one primal field, the size of this vector is one,otherwise the size is the
   *number of primal fields in the problem.
   */
  std::vector<FESystem<dim> *> FESet;
  /*A vector of all the constraint sets in the problem. A constraint set is a
   *map which holds the mapping between the degrees of freedom and the
   *corresponding degree of freedom constraints. Currently the type of
   *constraints stored are either Dirichlet boundary conditions or hanging node
   *constraints for adaptive meshes.
   */
  std::vector<const AffineConstraints<double> *> constraintsDirichletSet,
    constraintsOtherSet;
  /*A vector of all the degree of freedom objects is the problem. A degree of
   *freedom object handles the serial/parallel distribution of the degrees of
   *freedom for all the primal fields in the problem.*/
  std::vector<const DoFHandler<dim> *> dofHandlersSet;

  /*A vector of the locally relevant degrees of freedom. Locally relevant degrees of
   *freedom in a parallel implementation is a collection of the degrees of freedom owned
   *by the current processor and the surrounding ghost nodes which are required for the
   *field computations in this processor.
   */
  std::vector<const IndexSet *> locally_relevant_dofsSet;
  /*Copies of constraintSet elements, but stored as non-const to enable application of
   * constraints.*/
  std::vector<AffineConstraints<double> *> constraintsDirichletSet_nonconst,
    constraintsOtherSet_nonconst;
  /*Copies of dofHandlerSet elements, but stored as non-const.*/
  std::vector<DoFHandler<dim> *> dofHandlersSet_nonconst;
  /*Copies of locally_relevant_dofsSet elements, but stored as non-const.*/
  std::vector<IndexSet *> locally_relevant_dofsSet_nonconst;
  /*Vector all the solution vectors in the problem. In a multi-field problem, each primal
   * field has a solution vector associated with it.*/
  std::vector<vectorType *> solutionSet;
  /*Vector all the residual (RHS) vectors in the problem. In a multi-field problem, each
   * primal field has a residual vector associated with it.*/
  std::vector<vectorType *> residualSet;
  /*Vector of parallel solution transfer objects. This is used only when adaptive meshing
   * is enabled.*/
  std::vector<parallel::distributed::SolutionTransfer<dim, vectorType> *> soltransSet;

  // matrix free objects
  /*Object of class MatrixFree<dim>. This is primarily responsible for all the
   *base matrix free functionality of this MatrixFreePDE<dim> class. Refer to
   *deal.ii documentation of MatrixFree<dim> class for details.
   */
  MatrixFree<dim, double> matrixFreeObject;
  /*Vector to store the inverse of the mass matrix diagonal for scalar fields.
   * Due to the choice of spectral elements with Guass-Lobatto quadrature, the
   * mass matrix is diagonal.*/
  vectorType invMscalar;
  /*Vector to store the inverse of the mass matrix diagonal for vector fields.
   * Due to the choice of spectral elements with Guass-Lobatto quadrature, the
   * mass matrix is diagonal.*/
  vectorType invMvector;
  /*Vector to store the solution increment. This is a temporary vector used
   * during implicit solves of the Elliptic fields.*/
  vectorType dU_vector, dU_scalar;

  // matrix free methods
  /*Current field index*/
  unsigned int currentFieldIndex;
  /*Method to compute the inverse of the mass matrix*/
  void
  computeInvM();

  /*Method to compute an explicit timestep*/
  void
  updateExplicitSolution(unsigned int fieldIndex);

  /*Method to compute an implicit timestep*/
  bool
  updateImplicitSolution(unsigned int fieldIndex, unsigned int nonlinear_it_index);

  /*Method to apply boundary conditions*/
  void
  applyBCs(unsigned int fieldIndex);

  /**
   * \brief Compute element volume for the triangulation
   */
  void
  compute_element_volume();

  /**
   * \brief Vector that stores element volumes
   */
  dealii::AlignedVector<dealii::VectorizedArray<double>> element_volume;

  /*Method to compute the right hand side (RHS) residual vectors*/
  void
  computeExplicitRHS();
  void
  computeNonexplicitRHS();

  // virtual methods to be implemented in the derived class
  /*Method to calculate LHS(implicit solve)*/
  void
  getLHS(const MatrixFree<dim, double>               &data,
         vectorType                                  &dst,
         const vectorType                            &src,
         const std::pair<unsigned int, unsigned int> &cell_range) const;

  bool generatingInitialGuess;
  void
  getLaplaceLHS(const MatrixFree<dim, double>               &data,
                vectorType                                  &dst,
                const vectorType                            &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const;

  void
  setNonlinearEqInitialGuess();
  void
  computeLaplaceRHS(unsigned int fieldIndex);
  void
  getLaplaceRHS(const MatrixFree<dim, double>               &data,
                vectorType                                  &dst,
                const vectorType                            &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const;

  /*Method to calculate RHS (implicit/explicit). This is an abstract method, so
   * every model which inherits MatrixFreePDE<dim> has to implement this
   * method.*/
  void
  getExplicitRHS(const MatrixFree<dim, double>               &data,
                 std::vector<vectorType *>                   &dst,
                 const std::vector<vectorType *>             &src,
                 const std::pair<unsigned int, unsigned int> &cell_range) const;

  void
  getNonexplicitRHS(const MatrixFree<dim, double>               &data,
                    std::vector<vectorType *>                   &dst,
                    const std::vector<vectorType *>             &src,
                    const std::pair<unsigned int, unsigned int> &cell_range) const;

  virtual void
  explicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double>             element_volume) const = 0;

  virtual void
  nonExplicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double>             element_volume) const = 0;

  virtual void
  equationLHS([[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                                        &variable_list,
              [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
              [[maybe_unused]] const VectorizedArray<double> element_volume) const = 0;

  virtual void
  postProcessedFields(
    [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
      &variable_list,
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &pp_variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double>             element_volume) const {};
  void
  computePostProcessedFields(std::vector<vectorType *> &postProcessedSet);

  void
  getPostProcessedFields(const MatrixFree<dim, double>               &data,
                         std::vector<vectorType *>                   &dst,
                         const std::vector<vectorType *>             &src,
                         const std::pair<unsigned int, unsigned int> &cell_range);

  // methods to apply dirichlet BC's
  /*Map of degrees of freedom to the corresponding Dirichlet boundary
   * conditions, if any.*/
  std::vector<std::map<types::global_dof_index, double> *> valuesDirichletSet;
  /*Virtual method to mark the boundaries for applying Dirichlet boundary
   * conditions.  This is usually expected to be provided by the user.*/
  void
  markBoundaries(parallel::distributed::Triangulation<dim> &) const;
  /** Method for applying Dirichlet boundary conditions.*/
  void
  applyDirichletBCs();

  /** Method for applying Neumann boundary conditions.*/
  void
  applyNeumannBCs();

  // Methods to apply periodic BCs
  void
  setPeriodicity();
  void
  setPeriodicityConstraints(AffineConstraints<double> *, const DoFHandler<dim> *) const;

  /**
   * \brief Set constraints to pin the solution to 0 at a certain vertex. This is
   * automatically done at the origin if no value terms are detected in your dependencies
   * in a time_independent or implicit solve.
   *
   * \param constraints The constraint set.
   * \param dof_handler The list of the degrees of freedom.
   * \param target_point The point where the solution is constrained. This is the origin
   * by default.
   */
  void
  set_rigid_body_mode_constraints(AffineConstraints<double> *constraints,
                                  const DoFHandler<dim>     *dof_handler,
                                  const Point<dim> target_point = Point<dim>()) const;

  // methods to apply initial conditions
  /*Virtual method to apply initial conditions.  This is usually expected to be
   * provided by the user in IBVP (Initial Boundary Value Problems).*/

  void
  applyInitialConditions();

  // --------------------------------------------------------------------------
  // Methods for saving and loading checkpoints
  // --------------------------------------------------------------------------

  void
  save_checkpoint();

  void
  load_checkpoint_triangulation();
  void
  load_checkpoint_fields();
  void
  load_checkpoint_time_info();

  void
  move_file(const std::string &, const std::string &);

  void
  verify_checkpoint_file_exists(const std::string &filename);

  // --------------------------------------------------------------------------
  // Nucleation methods and variables
  // --------------------------------------------------------------------------
  // Vector of all the nuclei seeded in the problem
  std::vector<nucleus<dim>> nuclei;

  // Method to get a list of new nuclei to be seeded
  void
  updateNucleiList();
  std::vector<nucleus<dim>>
  getNewNuclei();
  void
  getLocalNucleiList(std::vector<nucleus<dim>> &newnuclei) const;
  void
  safetyCheckNewNuclei(std::vector<nucleus<dim>>  newnuclei,
                       std::vector<unsigned int> &conflict_ids);
  void
  refineMeshNearNuclei(std::vector<nucleus<dim>> newnuclei);
  double
  weightedDistanceFromNucleusCenter(const Point<dim, double>   center,
                                    const std::vector<double> &semiaxes,
                                    const Point<dim, double>   q_point_loc,
                                    const unsigned int         var_index) const;
  VectorizedArray<double>
  weightedDistanceFromNucleusCenter(const Point<dim, double>                  center,
                                    const std::vector<double>                &semiaxes,
                                    const Point<dim, VectorizedArray<double>> q_point_loc,
                                    const unsigned int var_index) const;

  // Method to obtain the nucleation probability for an element, nontrival case
  // must be implemented in the subsclass
  [[nodiscard]] virtual double
  getNucleationProbability(variableValueContainer,
                           double,
                           Point<dim>,
                           [[maybe_unused]] unsigned int variable_index) const
  {
    return 0.0;
  };

  // utility functions
  /*Returns index of given field name if exists, else throw error.*/
  unsigned int
  getFieldIndex(std::string _name);

  /*Method to compute the integral of a field.*/
  void
  computeIntegral(double                   &integratedField,
                  int                       index,
                  std::vector<vectorType *> variableSet);

  // variables for time dependent problems
  /*Flag used to see if invM, time stepping in run(), etc are necessary*/
  bool isTimeDependentBVP;
  /*Flag used to mark problems with Elliptic fields.*/
  bool isEllipticBVP;

  bool hasExplicitEquation;
  bool hasNonExplicitEquation;
  //
  double       currentTime;
  unsigned int currentIncrement, currentOutput, currentCheckpoint,
    current_grain_reassignment;

  /*Timer and logging object*/
  mutable TimerOutput computing_timer;

  bool first_integrated_var_output_complete;

  // Methods and variables for integration
  double       integrated_var;
  unsigned int integral_index;
  std::mutex   assembler_lock;

  /*AMR methods*/
  AdaptiveRefinement<dim, degree> AMR;
};

template class MatrixFreePDE<2, 1>;
template class MatrixFreePDE<3, 1>;

template class MatrixFreePDE<2, 2>;
template class MatrixFreePDE<3, 2>;

template class MatrixFreePDE<3, 3>;
template class MatrixFreePDE<2, 3>;

template class MatrixFreePDE<3, 4>;
template class MatrixFreePDE<2, 4>;

template class MatrixFreePDE<3, 5>;
template class MatrixFreePDE<2, 5>;

template class MatrixFreePDE<3, 6>;
template class MatrixFreePDE<2, 6>;

#endif
