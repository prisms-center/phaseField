// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/user_inputs/miscellaneous_parameters.h>
#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
MiscellaneousParameters::declare(dealii::ParameterHandler &parameter_handler,
                                 unsigned int              n_subsections)
{
  parameter_handler.enter_subsection("miscellaneous");
  {
    declare_entry(parameter_handler,
                  "random seed",
                  "2025",
                  Patterns::UnsignedInteger(),
                  "The random seed for the simulation. This is used to initialize the "
                  "random number generator.",
                  {},
                  "Our RNG engine adds the MPI process to the random seed. This means "
                  "that you must run a simulation with the same number of MPI processes "
                  "to get the same result when using RNG.");
  }
  parameter_handler.leave_subsection();
};

void
MiscellaneousParameters::assign(dealii::ParameterHandler &parameter_handler,
                                unsigned int              n_subsections)
{
  parameter_handler.enter_subsection("miscellaneous");
  {
    set_random_seed((unsigned int) (parameter_handler.get_integer("random seed")));
  }
  parameter_handler.leave_subsection();
};

void
MiscellaneousParameters::validate(const std::vector<FieldAttributes> &field_attributes,
                                  const std::vector<SolveBlock> &solve_blocks) const {
  // TODO: Do this later
};

void
MiscellaneousParameters::set_random_seed(const unsigned int &_random_seed)
{
  random_seed = _random_seed;
  rng.seed(random_seed + dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
}

PRISMS_PF_END_NAMESPACE
