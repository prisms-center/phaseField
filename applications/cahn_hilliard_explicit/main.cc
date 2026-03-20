// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/dependencies.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/phase_field_tools.h>
#include <prismspf/core/problem.h>
#include <prismspf/core/simulation_timer.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/input_file_reader.h>
#include <prismspf/user_inputs/miscellaneous_parameters.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <random>

using namespace prisms;

template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public PDEOperatorBase<dim, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using ScalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using VectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using VectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  /**
   * @brief Constructor.
   */
  explicit CustomPDE(const UserInputParameters<dim> &_user_inputs,
                     PhaseFieldTools<dim>           &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , McV(get_user_inputs().user_constants.get_model_constant_double("McV"))
    , KcV(get_user_inputs().user_constants.get_model_constant_double("KcV"))
    , WcV(get_user_inputs().user_constants.get_model_constant_double("WcV"))
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    const dealii::Tensor<1, dim> &mesh_size =
      get_user_inputs().spatial_discretization.rectangular_mesh.size;
    if (index == 0)
      {
        double center[12][3] = {
          {0.1, 0.3,  0},
          {0.8, 0.7,  0},
          {0.5, 0.2,  0},
          {0.4, 0.4,  0},
          {0.3, 0.9,  0},
          {0.8, 0.1,  0},
          {0.9, 0.5,  0},
          {0.0, 0.1,  0},
          {0.1, 0.6,  0},
          {0.5, 0.6,  0},
          {1,   1,    0},
          {0.7, 0.95, 0}
        };
        double rad[12] = {12, 14, 19, 16, 11, 12, 17, 15, 20, 10, 11, 14};
        double dist    = 0.0;
        double sdf     = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 12; i++)
          {
            dist = 0.0;
            for (unsigned int dir = 0; dir < dim; dir++)
              {
                double comp_diff = point[dir] - center[i][dir] * mesh_size[dir];
                dist += comp_diff * comp_diff;
              }
            dist = std::sqrt(dist) - rad[i];
            sdf  = std::min(sdf, dist);
          }
        scalar_value += 0.5 * (1.0 - std::tanh(sdf / 1.5));
      }
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_group_id) const override
  {
    if (solve_group_id == 0) // c
      {
        ScalarValue c = variable_list.template get_value<TensorRank::Scalar, OldOne>(0);
        ScalarGrad  mux =
          variable_list.template get_gradient<TensorRank::Scalar, OldOne>(1);

        ScalarValue eq_c  = c;
        ScalarGrad  eqx_c = -McV * sim_timer.get_timestep() * mux;

        variable_list.set_value_term(0, eq_c);
        variable_list.set_gradient_term(0, eqx_c);
      }
    else if (solve_group_id == 1) // mu
      {
        ScalarValue c = variable_list.template get_value<TensorRank::Scalar, Normal>(0);
        ScalarGrad  cx =
          variable_list.template get_gradient<TensorRank::Scalar, Normal>(0);

        ScalarValue fcV = WcV * c * (c - 1.0) * (c - 0.5);

        ScalarValue eq_mu  = fcV;
        ScalarGrad  eqx_mu = KcV * cx;

        variable_list.set_value_term(1, eq_mu);
        variable_list.set_gradient_term(1, eqx_mu);
      }
    else if (solve_group_id == 2) // pp
      {
        ScalarValue c = variable_list.template get_value<TensorRank::Scalar, Normal>(0);
        ScalarGrad  cx =
          variable_list.template get_gradient<TensorRank::Scalar, Normal>(0);

        ScalarValue f_tot  = 0.0;
        ScalarValue f_chem = c * c * c * c - 2.0 * c * c * c + c * c;
        ScalarValue f_grad = 0.0;

        for (unsigned int i = 0; i < dim; i++)
          {
            f_grad += 0.5 * KcV * cx[i] * cx[i];
          }
        f_tot = f_chem + f_grad;
        variable_list.set_value_term(2, f_tot);
      }
  }

  ScalarValue McV;
  ScalarValue KcV;
  ScalarValue WcV;
};

int
main(int argc, char *argv[])
{
  // Initialize MPI
  dealii::Utilities::MPI::MPI_InitFinalize
    mpi_init(argc, argv, dealii::numbers::invalid_unsigned_int);

  // Restrict deal.II console printing
  dealii::deallog.depth_console(0);

  // Parse the command line options (if there are any) to get the name of the input
  // file
  ParseCMDOptions cli_options(argc, argv);
  std::string     parameters_filename = cli_options.get_parameters_filename();

  constexpr unsigned int dim    = 2; // TODO change to 3 (original app)
  constexpr unsigned int degree = 2; // TODO change to 1 (original app)

  std::vector<FieldAttributes> fields = {FieldAttributes("c"),
                                         FieldAttributes("mu"),
                                         FieldAttributes("f_tot")};
  std::vector<SolveGroup>      solve_groups;
  SolveGroup                   c_group;
  c_group.id               = 0;
  c_group.solve_type       = Explicit;
  c_group.solve_timing     = Initialized;
  c_group.field_indices    = {0};
  c_group.dependencies_rhs = make_dependency_set(fields, {"old_1(c)", "grad(old_1(mu))"});

  SolveGroup mu_group;
  mu_group.id               = 1;
  mu_group.solve_type       = Explicit;
  mu_group.solve_timing     = Uninitialized;
  mu_group.field_indices    = {1};
  mu_group.dependencies_rhs = make_dependency_set(fields, {"c", "grad(c)"});

  SolveGroup pp_group;
  pp_group.id               = 2;
  pp_group.solve_type       = Explicit;
  pp_group.solve_timing     = PostProcess;
  pp_group.field_indices    = {2};
  pp_group.dependencies_rhs = make_dependency_set(fields, {"c", "grad(c)"});

  solve_groups.push_back(c_group);
  solve_groups.push_back(mu_group);
  solve_groups.push_back(pp_group);

  UserInputParameters<dim>       user_inputs(parameters_filename);
  PhaseFieldTools<dim>           pf_tools;
  CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);
  Problem<dim, degree, double>   problem(fields,
                                       solve_groups,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);
  problem.solve();

  return 0;
}
