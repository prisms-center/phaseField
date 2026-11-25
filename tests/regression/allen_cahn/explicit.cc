#include <prismspf/core/problem.h>

#include <prismspf/utilities/utilities.h>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

// TODO: Move to general purpose header
template <unsigned int dim, unsigned int degree, typename RealType>
struct RegressionTestConfig
{
  using real_type                            = RealType;
  static constexpr unsigned int dim_value    = dim;
  static constexpr unsigned int degree_value = degree;
};

using RegressionTestTypes = std::tuple<RegressionTestConfig<2, 2, float>>;

TEMPLATE_LIST_TEST_CASE("Explicit Allen-Cahn Regression",
                        "[regression][explicit][allen_cahn]",
                        RegressionTestTypes)
{
  using namespace prismspf;

  // Grab the template parameters
  using number = typename TestType::real_type;

  constexpr unsigned int dim    = TestType::dim_value;
  constexpr unsigned int degree = TestType::degree_value;

  // TODO: How and where do I initialize MPI

  // Define the PDE
  class ExplicitAllenCahn : public PDEOperatorBase<dim, degree, number>
  {
  public:
    ExplicitAllenCahn(const UserInputParameters<dim> &_user_inputs,
                      const PhaseFieldTools<dim>     &_pf_tools)
      : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools) {};

    void
    set_initial_condition([[maybe_unused]] const unsigned int       &index,
                          [[maybe_unused]] const unsigned int       &component,
                          [[maybe_unused]] const dealii::Point<dim> &point,
                          [[maybe_unused]] number                   &scalar_value,
                          [[maybe_unused]] number &vector_component_value) const override
    {}

    void
    compute_rhs(FieldContainer<dim, degree, number> &variable_list,
                const SimulationTimer               &sim_timer,
                unsigned int                         solve_block_id) const override
    {}
  };

  // Define the fields and solve blocks
  std::vector<FieldAttributes> fields = {FieldAttributes("phi")};

  SolveBlock phi_block;
  phi_block.id            = 0;
  phi_block.solve_type    = Explicit;
  phi_block.solve_timing  = Initialized;
  phi_block.field_indices = {0};
  phi_block.dependencies_rhs =
    make_dependency_set(fields, {"old_1(phi)", "grad(old_1(phi))"});

  std::vector<SolveBlock> solve_blocks({phi_block});

  // Define the parameters
  UserInputParameters<dim> user_inputs;
  user_inputs.spatial_discretization.global_refinement = 4;
  user_inputs.spatial_discretization.rectangular_mesh.size =
    make_tensor<dim, double>({10.0, 10.0, 10.0});

  // Other miscellaneous things
  PhaseFieldTools<dim>         pf_tools;
  ExplicitAllenCahn            pde_operator(user_inputs, pf_tools);
  Problem<dim, degree, number> problem(fields,
                                       solve_blocks,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);

  // Init system
  problem.init_system();

  // Solve a few steps and compare vector outputs
  REQUIRE(problem.get_solve_context().get_dof_manager().get_total_dofs() == 1);
}
