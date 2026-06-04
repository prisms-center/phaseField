// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_selector.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <map>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that stores relevant linear solve information of a certain solve block.
 */
struct LinearSolverParameters
{
  // Solver type. richardson|cg|bicgstab|gmres|fgmres|minres
  std::string solver_type = "cg";

  // Solver tolerance
  double tolerance = Defaults::tolerance;

  // Solver tolerance type
  SolverToleranceType tolerance_type = SolverToleranceType::RMSEPerField;

  // Max number of iterations for the linear solve
  unsigned int max_iterations = Defaults::iterations;

  // Preconditioner
  PreconditionerType preconditioner = PreconditionerType::None;

  dealii::PreconditionChebyshev<>::AdditionalData chebyshev_parameters;

  dealii::SolverRichardson<>::AdditionalData richardson_parameters;
  dealii::SolverBicgstab<>::AdditionalData   bicgstab_parameters;
  dealii::SolverGMRES<>::AdditionalData      gmres_parameters;

  // MinRes and CG do not have additional parameters, so we do not need to store them here
  // dealii::SolverMinRes<>::AdditionalData minres_parameters;
  // dealii::SolverCG<>::AdditionalData     cg_parameters;
  // FGMRES additional parameters are a subset of GMRES parameters, so we can just use
  // gmres_parameters for both solvers
  // dealii::SolverFGMRES<>::AdditionalData     fgmres_parameters;

  // The minimum multigrid level
  unsigned int min_mg_level = 0;
};

/**
 * @brief Struct that holds linear solver parameters.
 */
struct LinearSolveParameters
{
  /**
   * @brief Postprocess and validate parameters.
   */
  void
  validate();

  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler,
                     unsigned int              max_criteria = 5) const;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler,
                    unsigned int              max_criteria = 5);

  // Map of linear solve parameters for fields that require them
  std::map<unsigned int, LinearSolverParameters> linear_solvers;
};

inline void
LinearSolveParameters::validate()
{
  // Nothing to do here for now
}

inline void
LinearSolveParameters::print_parameter_summary() const
{
  if (!linear_solvers.empty())
    {
      ConditionalOStreams::pout_summary()
        << "================================================\n"
        << "  Linear Solve Parameters\n"
        << "================================================\n";

      for (const auto &[index, linear_solver_parameters] : linear_solvers)
        {
          ConditionalOStreams::pout_summary()
            << "Index: " << index << "\n"
            << "  Tolerance: " << linear_solver_parameters.tolerance << "\n"
            << "  Type: " << to_string(linear_solver_parameters.tolerance_type) << "\n"
            << "  Max iterations: " << linear_solver_parameters.max_iterations << "\n"
            << "  Preconditioner: " << to_string(linear_solver_parameters.preconditioner)
            << "\n";

          if (linear_solver_parameters.preconditioner == PreconditionerType::GMG)
            {
              // ConditionalOStreams::pout_summary();
            }
        }

      ConditionalOStreams::pout_summary() << "\n" << std::flush;
    }
}

inline void
LinearSolveParameters::declare_parameters(dealii::ParameterHandler &parameter_handler,
                                          unsigned int              max_criteria) const
{
  for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
    {
      // For linear solves
      std::string subsection_text =
        "linear solver parameters: " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        parameter_handler.declare_entry(
          "solver_ids",
          "",
          dealii::Patterns::Anything(),
          "The ids of the solvers that will use these settings.");
        parameter_handler.declare_entry(
          "solver type",
          "cg",
          dealii::Patterns::Selection(dealii::SolverSelector<>::get_solver_names()),
          "The type of iterative solver to use for linear solves.");
        parameter_handler.declare_entry(
          "tolerance type",
          "AbsoluteResidual",
          dealii::Patterns::Selection("AbsoluteResidual|RMSEPerField|IntegratedPerField|"
                                      "RMSETotal|IntegratedTotal"),
          "The tolerance type for the linear solver.");
        parameter_handler.declare_entry("tolerance value",
                                        "1.0e-10",
                                        dealii::Patterns::Double(DBL_MIN, DBL_MAX),
                                        "The value of for the linear solver tolerance.");
        parameter_handler.declare_entry(
          "max iterations",
          "100",
          dealii::Patterns::Integer(1, INT_MAX),
          "The maximum number of linear solver iterations before the loop "
          "is stopped.");
        parameter_handler.declare_entry(
          "preconditioner type",
          "None",
          dealii::Patterns::Selection(
            "|None|Chebyshev|GMG|none|chebyshev|gmg|MG|mg|multigrid"),
          "The preconditioner type for the linear solver.");
        parameter_handler.enter_subsection("Chebyshev");
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
        parameter_handler.leave_subsection();

        parameter_handler.enter_subsection("Richardson");
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
        parameter_handler.leave_subsection();
        parameter_handler.enter_subsection("BiCGStab");
        {
          parameter_handler.declare_entry("exact residual",
                                          "true",
                                          dealii::Patterns::Bool(),
                                          "Flag for exact computation of residual.");
          parameter_handler.declare_entry("breakdown",
                                          std::to_string(
                                            std::numeric_limits<double>::min()),
                                          dealii::Patterns::Double(),
                                          "Breakdown threshold.");
        }
        parameter_handler.leave_subsection();
        parameter_handler.enter_subsection("GMRES");
        {
          parameter_handler.declare_entry(
            "max basis size",
            "30",
            dealii::Patterns::Integer(1, INT_MAX),
            "The maximum size of the Krylov basis used in GMRES before restarting.");
          parameter_handler.declare_entry(
            "orthogonalization strategy",
            "delayed_classical_gram_schmidt",
            dealii::Patterns::Selection("classical_gram_schmidt|modified_gram_schmidt|"
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
        parameter_handler.leave_subsection();

        parameter_handler.declare_entry("min mg level",
                                        "0",
                                        dealii::Patterns::Integer(0, INT_MAX),
                                        "The minimum multigrid level.");
        parameter_handler.declare_alias("tolerance value", "tolerance");
        parameter_handler.declare_alias("solver_ids", "solve blocks");
        parameter_handler.declare_alias("solver_ids", "solve_blocks");
        parameter_handler.declare_alias("solver_ids", "solve block ids");
        parameter_handler.declare_alias("solver_ids", "solve_block_ids");
        parameter_handler.declare_alias("solver type", "solver_type");
        parameter_handler.declare_alias("solver type", "linear solver type");
        parameter_handler.declare_alias("solver type", "linear_solver_type");
        parameter_handler.declare_alias("solver type", "linear solver");
        parameter_handler.declare_alias("solver type", "linear_solver");
        parameter_handler.declare_alias("solver type", "type");
        parameter_handler.declare_alias("preconditioner type", "preconditioner_type");
        parameter_handler.declare_alias("preconditioner type", "preconditioner");
      }
      parameter_handler.leave_subsection();
    }
}

inline void
LinearSolveParameters::assign_parameters(dealii::ParameterHandler &parameter_handler,
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

        LinearSolverParameters linear_solver_parameters;

        // Set the linear solver type
        linear_solver_parameters.solver_type = parameter_handler.get("solver type");

        // Set the tolerance type
        static const std::map<std::string, SolverToleranceType> tolerance_types = {
          {"AbsoluteResidual",   AbsoluteResidual  },
          {"RMSEPerField",       RMSEPerField      },
          {"IntegratedPerField", IntegratedPerField},
          {"RMSETotal",          RMSETotal         },
          {"IntegratedTotal",    IntegratedTotal   },
        };
        const std::string type_string = parameter_handler.get("tolerance type");
        const auto       &iter        = tolerance_types.find(type_string);
        if (iter != tolerance_types.end())
          {
            linear_solver_parameters.tolerance_type = iter->second;
          }
        else
          {
            AssertThrow(false,
                        dealii::ExcMessage("Invalid tolerance type: " + type_string));
          }

        // Set the tolerance value
        linear_solver_parameters.tolerance =
          parameter_handler.get_double("tolerance value");

        // Set the maximum number of iterations
        linear_solver_parameters.max_iterations =
          static_cast<unsigned int>(parameter_handler.get_integer("max iterations"));

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
        linear_solver_parameters.preconditioner =
          preconditioner_map.at(parameter_handler.get("preconditioner type"));

        parameter_handler.enter_subsection("Chebyshev");
        {
          linear_solver_parameters.chebyshev_parameters.degree =
            static_cast<unsigned int>(parameter_handler.get_integer("smoother degree"));

          linear_solver_parameters.chebyshev_parameters.smoothing_range =
            parameter_handler.get_double("smoothing range");

          linear_solver_parameters.chebyshev_parameters.eig_cg_n_iterations =
            static_cast<unsigned int>(
              parameter_handler.get_integer("eigenvalue cg iterations"));
        }
        parameter_handler.leave_subsection();
        parameter_handler.enter_subsection("Richardson");
        {
          linear_solver_parameters.richardson_parameters.omega =
            parameter_handler.get_double("omega");
          linear_solver_parameters.richardson_parameters.use_preconditioned_residual =
            parameter_handler.get_bool("use preconditioned residual");
        }
        parameter_handler.leave_subsection();
        parameter_handler.enter_subsection("BiCGStab");
        {
          linear_solver_parameters.bicgstab_parameters.exact_residual =
            parameter_handler.get_bool("exact residual");
          linear_solver_parameters.bicgstab_parameters.breakdown =
            parameter_handler.get_double("breakdown");
        }
        parameter_handler.leave_subsection();
        parameter_handler.enter_subsection("GMRES");
        {
          using OStrat = dealii::LinearAlgebra::OrthogonalizationStrategy;
          static const std::map<std::string, OStrat> orthogonalization_strategy_names = {
            {"classical_gram_schmidt",         OStrat::classical_gram_schmidt        },
            {"modified_gram_schmidt",          OStrat::modified_gram_schmidt         },
            {"delayed_classical_gram_schmidt", OStrat::delayed_classical_gram_schmidt}
          };
          linear_solver_parameters.gmres_parameters.orthogonalization_strategy =
            orthogonalization_strategy_names.at(
              parameter_handler.get("orthogonalization strategy"));
          linear_solver_parameters.gmres_parameters.right_preconditioning =
            parameter_handler.get_bool("right preconditioning");
          linear_solver_parameters.gmres_parameters.use_default_residual =
            parameter_handler.get_bool("use default residual");
          linear_solver_parameters.gmres_parameters.force_re_orthogonalization =
            parameter_handler.get_bool("force re-orthogonalization");
          linear_solver_parameters.gmres_parameters.batched_mode =
            parameter_handler.get_bool("batched mode");
        }
        parameter_handler.leave_subsection();

        linear_solver_parameters.min_mg_level =
          static_cast<unsigned int>(parameter_handler.get_integer("min mg level"));
        for (auto solver_id : solver_ids)
          {
            linear_solvers[static_cast<unsigned int>(solver_id)] =
              linear_solver_parameters;
          }
      }
      parameter_handler.leave_subsection();
    }
}

PRISMS_PF_END_NAMESPACE
