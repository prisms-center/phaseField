#include <prismspf/user_inputs/miscellaneous_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
MiscellaneousParameters::declare(dealii::ParameterHandler &parameter_handler,
                                 unsigned int              n_subsections)
{
  parameter_handler.enter_subsection("miscellaneous");
  {
    parameter_handler.declare_entry(
      "random seed",
      "2025",
      dealii::Patterns::Integer(0, INT_MAX),
      "The random seed for the simulation. "
      "This is used to initialize the random number generator.");
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
