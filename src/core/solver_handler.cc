// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/solver_handler.h>
#include <prismspf/core/timer.h>

#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <ostream>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
SolveBlock<dim, degree, number>::SolveBlock(
  const SolverContext<dim, degree, number> &_solver_context,
  Types::Index                              block_index)
  : concurrent_constant_solver(_solver_context, block_index)
  , concurrent_explicit_solver(_solver_context, block_index)
  , concurrent_explicit_postprocess_solver(_solver_context, block_index)
  , sequential_auxiliary_solver(_solver_context, block_index)
  , sequential_linear_solver(_solver_context, block_index)
  , sequential_self_nonlinear_solver(_solver_context, block_index)
  , sequential_co_nonlinear_solver(_solver_context, block_index)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
SolveBlock<dim, degree, number>::init()
{
  concurrent_constant_solver.init();

  concurrent_explicit_solver.init();

  concurrent_explicit_postprocess_solver.init();

  sequential_auxiliary_solver.init();

  sequential_linear_solver.init();

  sequential_self_nonlinear_solver.init();

  sequential_co_nonlinear_solver.init();
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolveBlock<dim, degree, number>::reinit()
{
  concurrent_constant_solver.reinit();

  concurrent_explicit_solver.reinit();

  concurrent_explicit_postprocess_solver.reinit();

  sequential_auxiliary_solver.reinit();

  sequential_linear_solver.reinit();

  sequential_self_nonlinear_solver.reinit();

  sequential_co_nonlinear_solver.reinit();
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolveBlock<dim, degree, number>::solve(unsigned int increment, bool update_postprocessed)
{
  if (increment == 0)
    {
      // Solve the auxiliary fields at the 0th step
      ConditionalOStreams::pout_base() << "  solving auxiliary variables...\n"
                                       << std::flush;
      Timer::start_section("Auxiliary solver");
      sequential_auxiliary_solver.solve();
      Timer::end_section("Auxiliary solver");

      // Solve the linear time-independent fields at the 0th step
      ConditionalOStreams::pout_base()
        << "  solving linear time-independent variables...\n"
        << std::flush;
      Timer::start_section("Nonexplicit linear solver");
      sequential_linear_solver.solve();
      Timer::end_section("Nonexplicit linear solver");

      // Solve the self-nonlinear time-independent fields at the 0th step
      ConditionalOStreams::pout_base()
        << "  solving self-nonlinear time-independent variables...\n"
        << std::flush;
      Timer::start_section("Nonexplicit self-nonlinear solver");
      sequential_self_nonlinear_solver.solve();
      Timer::end_section("Nonexplicit self-nonlinear solver");

      // Solve the co-nonlinear time-independent fields at the 0th step
      ConditionalOStreams::pout_base()
        << "  solving co-nonlinear time-independent variables...\n"
        << std::flush;
      Timer::start_section("Nonexplicit co-nonlinear solver");
      sequential_co_nonlinear_solver.solve();
      Timer::end_section("Nonexplicit co-nonlinear solver");

      // Solve the postprocessed fields at the 0th step
      ConditionalOStreams::pout_base() << "  solving postprocessed variables...\n"
                                       << std::flush;
      Timer::start_section("Postprocess solver");
      concurrent_explicit_postprocess_solver.solve();
      Timer::end_section("Postprocess solver");

      return;
    }

  Timer::start_section("Explicit solver");
  concurrent_explicit_solver.solve();
  Timer::end_section("Explicit solver");

  Timer::start_section("Nonexplicit auxiliary solver");
  sequential_auxiliary_solver.solve();
  Timer::end_section("Nonexplicit auxiliary solver");

  Timer::start_section("Nonexplicit linear solver");
  sequential_linear_solver.solve();
  Timer::end_section("Nonexplicit linear solver");

  Timer::start_section("Nonexplicit self-nonlinear solver");
  sequential_self_nonlinear_solver.solve();
  Timer::end_section("Nonexplicit self-nonlinear solver");

  Timer::start_section("Nonexplicit co-nonlinear solver");
  sequential_co_nonlinear_solver.solve();
  Timer::end_section("Nonexplicit co-nonlinear solver");

  if (update_postprocessed)
    {
      Timer::start_section("Postprocess solver");
      concurrent_explicit_postprocess_solver.solve();
      Timer::end_section("Postprocess solver");
    }
}

template <unsigned int dim, unsigned int degree, typename number>
SolverHandler<dim, degree, number>::SolverHandler(
  const SolverContext<dim, degree, number> &_solver_context)
{
  // Create a set of the solve blocks
  for (const auto &[var_index, variable] :
       _solver_context.get_user_inputs().get_variable_attributes())
    {
      solve_blocks.try_emplace(variable.get_solve_block(),
                               _solver_context,
                               variable.get_solve_block());
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverHandler<dim, degree, number>::init()
{
  Timer::start_section("Solver initialization");

  for (auto &[block_index, solve_block] : solve_blocks)
    {
      ConditionalOStreams::pout_base() << "Initializing solvers for solve block "
                                       << std::to_string(block_index) << " ...\n"
                                       << std::flush;
      solve_block.init();
    }

  Timer::end_section("Solver initialization");
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverHandler<dim, degree, number>::reinit()
{
  Timer::start_section("Solver reinitialization");

  for (auto &[block_index, solve_block] : solve_blocks)
    {
      solve_block.reinit();
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
      for (auto &[block_index, solve_block] : solve_blocks)
        {
          ConditionalOStreams::pout_base() << "solving 0th timestep for solve block "
                                           << std::to_string(block_index) << " ...\n"
                                           << std::flush;
          solve_block.solve(increment, update_postprocessed);
        }
      return;
    }

  for (auto &[block_index, solve_block] : solve_blocks)
    {
      solve_block.solve(increment, update_postprocessed);
    }
}

#include "core/solver_handler.inst"

PRISMS_PF_END_NAMESPACE
