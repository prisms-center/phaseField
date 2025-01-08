#ifndef explicit_solver_h
#define explicit_solver_h

/**
 * \brief This class handles the explicit solves of all explicit fields
 */
class explicitSolver
{
public:
  /**
   * \brief Constructor.
   */
  explicitSolver();

  /**
   * \brief Destructor.
   */
  ~explicitSolver();

  /**
   * \brief Initialize system.
   */
  void
  init_system();

  /**
   * \brief Reinitialize system.
   */
  void
  reinit_system();

  /**
   * \brief Solve a single update step.
   */
  void
  solve();
};

#endif