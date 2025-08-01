// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/solver_handler.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
SolverHandler<dim, degree, number>::SolverHandler(
  const SolverContext<dim, degree> &_solver_context)
{
  // Create a set of the solve blocks
  for (const auto &[index, variable] :
       _solver_context.get_user_inputs().get_variable_attributes())
    {
      solve_blocks.insert(variable.get_solve_block());
    }

  // Create the map of solvers according to each solve block
  for (const auto &solve_block : solve_blocks)
    {
      concurrent_constant_solver.try_emplace(solve_block, _solver_context, solve_block);
      concurrent_explicit_solver.try_emplace(solve_block, _solver_context, solve_block);
      concurrent_explicit_postprocess_solver.try_emplace(solve_block,
                                                         _solver_context,
                                                         solve_block);
      sequential_auxiliary_solver.try_emplace(solve_block, _solver_context, solve_block);
      sequential_linear_solver.try_emplace(solve_block, _solver_context, solve_block);
      sequential_self_nonlinear_solver.try_emplace(solve_block,
                                                   _solver_context,
                                                   solve_block);
      sequential_co_nonlinear_solver.try_emplace(solve_block,
                                                 _solver_context,
                                                 solve_block);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverHandler<dim, degree, number>::init()
{
  Timer::start_section("Solver initialization");

  for (const auto &solve_block : solve_blocks)
    {
      ConditionalOStreams::pout_base()
        << "initializing solvers for solve block " << std::to_string(solve_block)
        << " ...\n " << std::flush;

      ConditionalOStreams::pout_base()
        << "  trying to initialize concurrent constant solvers...\n"
        << std::flush;
      concurrent_constant_solver.at(solve_block).init();

      ConditionalOStreams::pout_base()
        << "  trying to initialize concurrent explicit solvers...\n"
        << std::flush;
      concurrent_explicit_solver.at(solve_block).init();

      ConditionalOStreams::pout_base()
        << "  trying to initialize concurrent explicit postprocess solvers...\n"
        << std::flush;
      concurrent_explicit_postprocess_solver.at(solve_block).init();

      ConditionalOStreams::pout_base()
        << "  trying to initialize sequential auxiliary solvers...\n"
        << std::flush;
      sequential_auxiliary_solver.at(solve_block).init();

      ConditionalOStreams::pout_base()
        << "  trying to initialize sequential linear solvers...\n"
        << std::flush;
      sequential_linear_solver.at(solve_block).init();

      ConditionalOStreams::pout_base()
        << "  trying to initialize sequential self-nonlinear solvers...\n"
        << std::flush;
      sequential_self_nonlinear_solver.at(solve_block).init();

      ConditionalOStreams::pout_base()
        << "  trying to initialize sequential co-nonlinear solvers...\n"
        << std::flush;
      sequential_co_nonlinear_solver.at(solve_block).init();
    }

  Timer::end_section("Solver initialization");
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverHandler<dim, degree, number>::reinit()
{
  Timer::start_section("Solver reinitialization");

  for (const auto &solve_block : solve_blocks)
    {
      ConditionalOStreams::pout_base()
        << "reinitializing solvers for solve block " << std::to_string(solve_block)
        << " ...\n " << std::flush;

      ConditionalOStreams::pout_base()
        << "  trying to reinitialize concurrent constant solvers...\n"
        << std::flush;
      concurrent_constant_solver.at(solve_block).reinit();

      ConditionalOStreams::pout_base()
        << "  trying to reinitialize concurrent explicit solvers...\n"
        << std::flush;
      concurrent_explicit_solver.at(solve_block).reinit();

      ConditionalOStreams::pout_base()
        << "  trying to reinitialize concurrent explicit postprocess solvers...\n"
        << std::flush;
      concurrent_explicit_postprocess_solver.at(solve_block).reinit();

      ConditionalOStreams::pout_base()
        << "  trying to reinitialize sequential auxiliary solvers...\n"
        << std::flush;
      sequential_auxiliary_solver.at(solve_block).reinit();

      ConditionalOStreams::pout_base()
        << "  trying to reinitialize sequential linear solvers...\n"
        << std::flush;
      sequential_linear_solver.at(solve_block).reinit();

      ConditionalOStreams::pout_base()
        << "  trying to reinitialize sequential self-nonlinear solvers...\n"
        << std::flush;
      sequential_self_nonlinear_solver.at(solve_block).reinit();

      ConditionalOStreams::pout_base()
        << "  trying to reinitialize sequential co-nonlinear solvers...\n"
        << std::flush;
      sequential_co_nonlinear_solver.at(solve_block).reinit();
    }

  Timer::end_section("Solver reinitialization");
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverHandler<dim, degree, number>::solve(unsigned int increment,
                                          bool         update_postprocessed)
{
  if (increment == 0)
    {
      for (const auto &solve_block : solve_blocks)
        {
          ConditionalOStreams::pout_base()
            << "solving 0th timestep for solve block " << std::to_string(solve_block)
            << " ...\n " << std::flush;

          // Solve the auxiliary fields at the 0th step
          ConditionalOStreams::pout_base() << "  solving auxiliary variables...\n"
                                           << std::flush;
          Timer::start_section("Auxiliary solver");
          sequential_auxiliary_solver.at(solve_block).solve();
          Timer::end_section("Auxiliary solver");

          // Solve the linear time-independent fields at the 0th step
          ConditionalOStreams::pout_base()
            << "  solving linear time-independent variables...\n"
            << std::flush;
          Timer::start_section("Nonexplicit linear solver");
          sequential_linear_solver.at(solve_block).solve();
          Timer::end_section("Nonexplicit linear solver");

          // Solve the self-nonlinear time-independent fields at the 0th step
          ConditionalOStreams::pout_base()
            << "  solving self-nonlinear time-independent variables...\n"
            << std::flush;
          Timer::start_section("Nonexplicit self-nonlinear solver");
          sequential_self_nonlinear_solver.at(solve_block).solve();
          Timer::end_section("Nonexplicit self-nonlinear solver");

          // Solve the co-nonlinear time-independent fields at the 0th step
          ConditionalOStreams::pout_base()
            << "  solving co-nonlinear time-independent variables...\n"
            << std::flush;
          Timer::start_section("Nonexplicit co-nonlinear solver");
          sequential_co_nonlinear_solver.at(solve_block).solve();
          Timer::end_section("Nonexplicit co-nonlinear solver");

          if (update_postprocessed)
            {
              // Solve the postprocessed fields at the 0th step
              ConditionalOStreams::pout_base() << "  solving postprocessed variables...\n"
                                               << std::flush;
              Timer::start_section("Postprocess solver");
              concurrent_explicit_postprocess_solver.at(solve_block).solve();
              Timer::end_section("Postprocess solver");
            }
        }
      return;
    }

  for (const auto &solve_block : solve_blocks)
    {
      Timer::start_section("Explicit solver");
      concurrent_explicit_solver.at(solve_block).solve();
      Timer::end_section("Explicit solver");

      Timer::start_section("Nonexplicit auxiliary solver");
      sequential_auxiliary_solver.at(solve_block).solve();
      Timer::end_section("Nonexplicit auxiliary solver");

      Timer::start_section("Nonexplicit linear solver");
      sequential_linear_solver.at(solve_block).solve();
      Timer::end_section("Nonexplicit linear solver");

      Timer::start_section("Nonexplicit self-nonlinear solver");
      sequential_self_nonlinear_solver.at(solve_block).solve();
      Timer::end_section("Nonexplicit self-nonlinear solver");

      Timer::start_section("Nonexplicit co-nonlinear solver");
      sequential_co_nonlinear_solver.at(solve_block).solve();
      Timer::end_section("Nonexplicit co-nonlinear solver");

      if (update_postprocessed)
        {
          Timer::start_section("Postprocess solver");
          concurrent_explicit_postprocess_solver.at(solve_block).solve();
          Timer::end_section("Postprocess solver");
        }
    };
}

#include "core/solver_handler.inst"

PRISMS_PF_END_NAMESPACE