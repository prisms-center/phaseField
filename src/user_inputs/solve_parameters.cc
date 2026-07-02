#include <prismspf/user_inputs/solve_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
LinearSolverParameters::declare(dealii::ParameterHandler &parameter_handler,
                                unsigned int              max_criteria) const
{
  parameter_handler.declare_entry(
    "solver type",
    "cg",
    dealii::Patterns::Selection(dealii::SolverSelector<>::get_solver_names()),
    "The type of iterative solver to use for linear solves.");
  declare_aliases(parameter_handler,
                  "solver type",
                  std::vector {"solver_type",
                               "linear solver type",
                               "linear_solver_type",
                               "linear solver",
                               "linear_solver",
                               "type"});

  parameter_handler.declare_entry("tolerance type",
                                  "AbsoluteResidual",
                                  dealii::Patterns::Selection(
                                    "AbsoluteResidual|RMSEPerField|IntegratedPerField|"
                                    "RMSETotal|IntegratedTotal"),
                                  "The tolerance type for the linear solver.");

  parameter_handler.declare_entry("tolerance value",
                                  "1.0e-10",
                                  dealii::Patterns::Double(0.0, DBL_MAX),
                                  "The value of for the linear solver tolerance.");
  parameter_handler.declare_alias("tolerance value", "tolerance");

  parameter_handler.declare_entry(
    "max iterations",
    "100",
    dealii::Patterns::Integer(1, INT_MAX),
    "The maximum number of linear solver iterations before the loop "
    "is stopped.");

  parameter_handler.declare_entry("mg depth",
                                  "1",
                                  dealii::Patterns::Integer(1, INT_MAX),
                                  "The depth of the multigrid hierarchy.");
  parameter_handler.declare_alias("mg depth", "mg_depth");

  parameter_handler.declare_entry(
    "preconditioner type",
    "None",
    dealii::Patterns::Selection("None|Chebyshev|GMG|none|chebyshev|gmg|MG|mg|multigrid"),
    "The preconditioner type for the linear solver.");
  declare_aliases(parameter_handler,
                  "preconditioner type",
                  std::vector {"preconditioner_type", "preconditioner"});

  // Now declare parameters for each of the solver's AdditionalData structures.
  parameter_handler.enter_subsection("Chebyshev");
  {
    declare_chebyshev(parameter_handler);
  }
  parameter_handler.leave_subsection();
  parameter_handler.enter_subsection("Richardson");
  {
    declare_richardson(parameter_handler);
  }
  parameter_handler.leave_subsection();
  parameter_handler.enter_subsection("BiCGStab");
  {
    declare_bicgstab(parameter_handler);
  }
  parameter_handler.leave_subsection();
  parameter_handler.enter_subsection("GMRES");
  {
    declare_gmres(parameter_handler);
  }
  parameter_handler.leave_subsection();
}

void
LinearSolverParameters::assign(dealii::ParameterHandler &parameter_handler,
                               unsigned int              max_criteria)
{
  // Set the linear solver type
  solver_type = parameter_handler.get("solver type");

  // Set the tolerance type
  static const std::map<std::string, SolverToleranceType> tolerance_types = {
    {"AbsoluteResidual",   AbsoluteResidual  },
    {"RMSEPerField",       RMSEPerField      },
    {"IntegratedPerField", IntegratedPerField},
    {"RMSETotal",          RMSETotal         },
    {"IntegratedTotal",    IntegratedTotal   },
  };
  const std::string type_string = parameter_handler.get("tolerance type");

  AssertThrow(tolerance_types.contains(type_string),
              dealii::ExcMessage("Invalid tolerance type: " + type_string));
  tolerance_type = tolerance_types.at(type_string);

  // Set the tolerance value
  tolerance = parameter_handler.get_double("tolerance value");

  // Set the maximum number of iterations
  max_iterations = (unsigned int) (parameter_handler.get_integer("max iterations"));

  mg_depth = (unsigned int) (parameter_handler.get_integer("mg depth"));

  // Set preconditioner type and related parameters
  static const std::map<std::string, PreconditionerType> preconditioner_map = {
    {"None",      None     },
    {"none",      None     },
    {"",          None     },
    {"Chebyshev", Chebyshev},
    {"chebyshev", Chebyshev},
    {"GMG",       GMG      },
    {"gmg",       GMG      },
    {"MG",        GMG      },
    {"mg",        GMG      },
    {"multigrid", GMG      }
  };
  preconditioner = preconditioner_map.at(parameter_handler.get("preconditioner type"));

  parameter_handler.enter_subsection("Chebyshev");
  {
    assign_chebyshev(parameter_handler);
  }
  parameter_handler.leave_subsection();
  parameter_handler.enter_subsection("Richardson");
  {
    assign_richardson(parameter_handler);
  }
  parameter_handler.leave_subsection();
  parameter_handler.enter_subsection("BiCGStab");
  {
    assign_bicgstab(parameter_handler);
  }
  parameter_handler.leave_subsection();
  parameter_handler.enter_subsection("GMRES");
  {
    assign_gmres(parameter_handler);
  }
  parameter_handler.leave_subsection();
}

void
LinearSolverParameters::validate(const std::vector<FieldAttributes> &field_attributes,
                                 const std::vector<SolveBlock>      &solve_blocks) const
{
  // TODO: Add validation here
}

void
LinearSolverParameters::declare_chebyshev(
  dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.declare_entry("smoothing range",
                                  "15.0",
                                  dealii::Patterns::Double(DBL_MIN, DBL_MAX),
                                  "The smoothing range for eigenvalues.");
  parameter_handler.declare_entry("smoother degree",
                                  "5",
                                  dealii::Patterns::Integer(1, INT_MAX),
                                  "The smoother polynomial degree.");
  parameter_handler.declare_entry(
    "eigenvalue cg iterations",
    "10",
    dealii::Patterns::Integer(1, INT_MAX),
    "The maximum number of CG iterations used to find the maximum eigenvalue.");
}

void
LinearSolverParameters::assign_chebyshev(dealii::ParameterHandler &parameter_handler)
{
  chebyshev_parameters.degree =
    (unsigned int) (parameter_handler.get_integer("smoother degree"));

  chebyshev_parameters.smoothing_range = parameter_handler.get_double("smoothing range");

  chebyshev_parameters.eig_cg_n_iterations =
    (unsigned int) (parameter_handler.get_integer("eigenvalue cg iterations"));
}

void
LinearSolverParameters::declare_richardson(
  dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.declare_entry("omega",
                                  "1.0",
                                  dealii::Patterns::Double(),
                                  "Damping factor.");
  parameter_handler.declare_entry("use preconditioned residual",
                                  "false",
                                  dealii::Patterns::Bool(),
                                  "Whether to use the preconditioned residual l2 "
                                  "norm in the stopping criterion.");
}

void
LinearSolverParameters::assign_richardson(dealii::ParameterHandler &parameter_handler)
{
  richardson_parameters.omega = parameter_handler.get_double("omega");
  richardson_parameters.use_preconditioned_residual =
    parameter_handler.get_bool("use preconditioned residual");
}

void
LinearSolverParameters::declare_bicgstab(
  dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.declare_entry("exact residual",
                                  "true",
                                  dealii::Patterns::Bool(),
                                  "Flag for exact computation of residual.");
  parameter_handler.declare_entry("breakdown",
                                  std::to_string(std::numeric_limits<double>::min()),
                                  dealii::Patterns::Double(),
                                  "Breakdown threshold.");
}

void
LinearSolverParameters::assign_bicgstab(dealii::ParameterHandler &parameter_handler)
{
  bicgstab_parameters.exact_residual = parameter_handler.get_bool("exact residual");
  bicgstab_parameters.breakdown      = parameter_handler.get_double("breakdown");
}

void
LinearSolverParameters::declare_gmres(dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.declare_entry(
    "max basis size",
    "30",
    dealii::Patterns::Integer(1, INT_MAX),
    "The maximum size of the Krylov basis used in GMRES before restarting.");
  parameter_handler.declare_entry("orthogonalization strategy",
                                  "delayed_classical_gram_schmidt",
                                  dealii::Patterns::Selection(
                                    "classical_gram_schmidt|modified_gram_schmidt|"
                                    "delayed_classical_gram_schmidt"),
                                  "The orthogonalization strategy to use in GMRES.");
  parameter_handler.declare_entry("right preconditioning",
                                  "false",
                                  dealii::Patterns::Bool(),
                                  "Whether to use right preconditioning.");
  parameter_handler.declare_entry(
    "use default residual",
    "true",
    dealii::Patterns::Bool(),
    "Whether to use the default residual computation in GMRES.");
  parameter_handler.declare_entry(
    "force re-orthogonalization",
    "false",
    dealii::Patterns::Bool(),
    "Whether to force re-orthogonalization of the Krylov basis in GMRES.");
  parameter_handler.declare_entry("batched mode",
                                  "false",
                                  dealii::Patterns::Bool(),
                                  "Whether to use batched mode in GMRES.");
}

void
LinearSolverParameters::assign_gmres(dealii::ParameterHandler &parameter_handler)
{
  using OStrat = dealii::LinearAlgebra::OrthogonalizationStrategy;
  static const std::map<std::string, OStrat> orthogonalization_strategy_names = {
    {"classical_gram_schmidt",         OStrat::classical_gram_schmidt        },
    {"modified_gram_schmidt",          OStrat::modified_gram_schmidt         },
    {"delayed_classical_gram_schmidt", OStrat::delayed_classical_gram_schmidt}
  };
  gmres_parameters.orthogonalization_strategy = orthogonalization_strategy_names.at(
    parameter_handler.get("orthogonalization strategy"));
  gmres_parameters.right_preconditioning =
    parameter_handler.get_bool("right preconditioning");
  gmres_parameters.use_default_residual =
    parameter_handler.get_bool("use default residual");
  gmres_parameters.force_re_orthogonalization =
    parameter_handler.get_bool("force re-orthogonalization");
  gmres_parameters.batched_mode = parameter_handler.get_bool("batched mode");
}

void
LinearSolveParameters::declare(dealii::ParameterHandler &parameter_handler,
                               unsigned int              max_criteria) const
{
  for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
    {
      std::string subsection_text =
        "linear solver parameters: " + std::to_string(criterion_id);

      parameter_handler.enter_subsection(subsection_text);
      {
        parameter_handler.declare_entry(
          "solver_ids",
          "",
          dealii::Patterns::List(dealii::Patterns::Anything(), 0, INT_MAX, ","),
          "The ids of the solvers that will use these settings.");
        declare_aliases(parameter_handler,
                        "solver_ids",
                        std::vector {"solve blocks",
                                     "solve_blocks",
                                     "solve block ids",
                                     "solve_block_ids",
                                     "solver ids"});

        LinearSolverParameters linear_solver;
        linear_solver.declare(parameter_handler);
      }
      parameter_handler.leave_subsection();
    }
}

void
LinearSolveParameters::assign(dealii::ParameterHandler &parameter_handler,
                              unsigned int              max_criteria)
{
  for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
    {
      // For linear solves
      std::string subsection_text =
        "linear solver parameters: " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        std::vector<int> solver_ids = dealii::Utilities::string_to_int(
          dealii::Utilities::split_string_list(parameter_handler.get("solver_ids")));

        LinearSolverParameters linear_solver;
        linear_solver.assign(parameter_handler);

        for (auto solver_id : solver_ids)
          {
            linear_solvers[(unsigned int) (solver_id)] = linear_solver;
          }
      }
      parameter_handler.leave_subsection();
    }
}

void
LinearSolveParameters::validate(const std::vector<FieldAttributes> &field_attributes,
                                const std::vector<SolveBlock>      &solve_blocks) const
{
  for (const auto &[solver_id, linear_solver] : linear_solvers)
    {
      linear_solver.validate(field_attributes, solve_blocks);
    }
}

void
NonlinearSolverParameters::declare(dealii::ParameterHandler &parameter_handler,
                                   unsigned int              max_criteria) const
{
  parameter_handler.declare_entry("max iterations",
                                  "100",
                                  dealii::Patterns::Integer(1, INT_MAX),
                                  "The maximum number of nonlinear solver "
                                  "iterations before the loop is stopped.");

  /*         parameter_handler.declare_entry(
        "tolerance type",
        "AbsoluteResidual",
        dealii::Patterns::Selection("AbsoluteResidual|RMSEPerField|IntegratedPerField|"
                                    "RMSETotal|IntegratedTotal"),
        "The tolerance type for the nonlinear solver."); */

  parameter_handler.declare_entry("tolerance value",
                                  "1.0e-10",
                                  dealii::Patterns::Double(0.0, DBL_MAX),
                                  "The value of for the nonlinear solver tolerance.");
  parameter_handler.declare_alias("tolerance value", "tolerance");

  /*         parameter_handler.declare_entry(
        "use backtracking line search",
        "true",
        dealii::Patterns::Bool(),
        "Whether to use a backtracking line-search to find the best "
        "choice of the damping coefficient.");
      parameter_handler.declare_entry(
        "step size modifier",
        "0.5",
        dealii::Patterns::Double(0.0, 1.0),
        "The constant that determines how much the step size decreases "
        "per backtrack. The 'tau' parameter.");
      parameter_handler.declare_entry(
        "residual decrease coefficient",
        "0.5",
        dealii::Patterns::Double(0.0, 1.0),
        "The constant that determines how much the residual must "
        "decrease to be accepted as sufficient. The 'c' parameter."); */

  parameter_handler.declare_entry(
    "step size",
    "1.0",
    dealii::Patterns::Double(0.0),
    "The constant damping value to be used if the backtrace "
    "line-search approach isn't used.");
}

void
NonlinearSolverParameters::assign(dealii::ParameterHandler &parameter_handler,
                                  unsigned int              max_criteria)
{
  max_iterations = (unsigned int) (parameter_handler.get_integer("max iterations"));

  step_length = parameter_handler.get_double("step size");

  tolerance_value = parameter_handler.get_double("tolerance value");
}

void
NonlinearSolverParameters::validate(const std::vector<FieldAttributes> &field_attributes,
                                    const std::vector<SolveBlock> &solve_blocks) const
{
  // TODO: Add validation here
}

void
NonlinearSolveParameters::declare(dealii::ParameterHandler &parameter_handler,
                                  unsigned int              max_criteria) const
{
  for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
    {
      std::string subsection_text =
        "newton solver parameters: " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        parameter_handler.declare_entry(
          "solver_ids",
          "",
          dealii::Patterns::List(dealii::Patterns::Anything(), 0, INT_MAX, ","),
          "The ids of the solvers that will use these settings.");
        declare_aliases(parameter_handler,
                        "solver_ids",
                        std::vector {"solve blocks",
                                     "solve_blocks",
                                     "solve block ids",
                                     "solve_block_ids",
                                     "solver ids"});

        NonlinearSolverParameters nonlinear_solver;
        nonlinear_solver.declare(parameter_handler);
      }
      parameter_handler.leave_subsection();
    }
}

void
NonlinearSolveParameters::assign(dealii::ParameterHandler &parameter_handler,
                                 unsigned int              max_criteria)
{
  for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
    {
      std::string subsection_text =
        "newton solver parameters: " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        std::vector<int> solver_ids = dealii::Utilities::string_to_int(
          dealii::Utilities::split_string_list(parameter_handler.get("solver_ids")));

        NonlinearSolverParameters nonlinear_solver;
        nonlinear_solver.assign(parameter_handler);

        for (auto solver_id : solver_ids)
          {
            nonlinear_solvers[(unsigned int) (solver_id)] = nonlinear_solver;
          }
      }
      parameter_handler.leave_subsection();
    }
}

void
NonlinearSolveParameters::validate(const std::vector<FieldAttributes> &field_attributes,
                                   const std::vector<SolveBlock>      &solve_blocks) const
{
  for (const auto &[solver_id, nonlinear_solver] : nonlinear_solvers)
    {
      nonlinear_solver.validate(field_attributes, solve_blocks);
    }
}

PRISMS_PF_END_NAMESPACE
