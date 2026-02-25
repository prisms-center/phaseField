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
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

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
    , m_well(1.0)
    , kappa(2.0)
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    number center[12][3] = {
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
    number rad[12] = {12, 14, 19, 16, 11, 12, 17, 15, 20, 10, 11, 14};
    number dist    = 0.0;
    for (unsigned int i = 0; i < 12; i++)
      {
        dist = 0.0;
        for (unsigned int dir = 0; dir < dim; dir++)
          {
            dist += (point[dir] - center[i][dir] * 100.0) *
                    (point[dir] - center[i][dir] * 100.0);
          }
        dist = std::sqrt(dist);

        scalar_value += 0.5 * (1.0 - std::tanh((dist - rad[i]) / 1.5));
      }
    scalar_value = std::min(scalar_value, static_cast<number>(1.0));
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_group_id) const override
  {
    ScalarValue n_val =
      variable_list.template get_value<TensorRank::Scalar, DependencyType::OldOne>(0);
    ScalarGrad n_grad =
      variable_list.template get_gradient<TensorRank::Scalar, DependencyType::OldOne>(0);

    if (solve_group_id == 1) // explicit
      {
        ScalarValue f_well = 4.0 * n_val * (n_val - 1.0) * (n_val - 0.5);
        ScalarValue eq_n   = n_val - sim_timer.get_timestep() * m_well * f_well;
        ScalarGrad  eqx_n  = -sim_timer.get_timestep() * kappa * m_well * n_grad;

        variable_list.set_value_term(0, eq_n);
        variable_list.set_gradient_term(0, eqx_n);
      }
    else if (solve_group_id == 2) // postprocess
      {
        ScalarValue f_tot = 0.0;
        ScalarValue f_chem =
          n_val * n_val * n_val * n_val - 2.0 * n_val * n_val * n_val + n_val * n_val;
        ScalarValue f_grad = 0.0;
        for (unsigned int i = 0; i < dim; i++)
          {
            f_grad += 0.5 * kappa * n_grad[i] * n_grad[i];
          }
        f_tot = f_chem + f_grad;

        variable_list.set_value_term(1,
                                     std::sqrt(n_grad[0] * n_grad[0] +
                                               n_grad[1] * n_grad[1]));
        variable_list.set_value_term(2, f_tot);
      }
  }

  number m_well;
  number kappa;
};

int
main(int argc, char *argv[])
{
  try
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

      constexpr unsigned int dim    = 3;
      constexpr unsigned int degree = 2;

      std::vector<FieldAttributes> field_attributes = {
        FieldAttributes("n", TensorRank::Scalar),
        FieldAttributes("mg_n", TensorRank::Scalar),
        FieldAttributes("f_tot", TensorRank::Scalar),
      };
      std::vector<SolveGroup> solve_groups;
      SolveGroup              exp_group(1,
                           Explicit,
                           Primary,
                                        {0},
                           make_dependency_set(field_attributes,
                                               std::set<std::string> {"n", "grad(n)"}));
      SolveGroup              pp_group(2,
                          Explicit,
                          PostProcess,
                                       {1, 2},
                          make_dependency_set(field_attributes,
                                              std::set<std::string> {"n", "grad(n)"}));
      solve_groups.push_back(exp_group);
      solve_groups.push_back(pp_group);

      InputFileReader prm(parameters_filename);
      const auto     &param_handler = prm.get_parameter_handler();
      param_handler.print_parameters("thing.prm", dealii::ParameterHandler::PRM);
      UserInputParameters<dim> user_inputs(parameters_filename);
      PhaseFieldTools<dim>     pf_tools;
      std::shared_ptr<PDEOperatorBase<dim, degree, double>> pde_operator =
        std::make_shared<CustomPDE<dim, degree, double>>(user_inputs, pf_tools);
      Problem<dim, degree, double> problem(field_attributes,
                                           solve_groups,
                                           user_inputs,
                                           pf_tools,
                                           pde_operator);
      problem.solve();
    }

  catch (std::exception &exc)
    {
      std::cerr << '\n'
                << '\n'
                << "----------------------------------------------------" << '\n';
      std::cerr << "Exception on: " << '\n'
                << exc.what() << '\n'
                << "Aborting!" << '\n'
                << "----------------------------------------------------" << '\n';
      return 1;
    }

  catch (...)
    {
      std::cerr << '\n'
                << '\n'
                << "----------------------------------------------------" << '\n';
      std::cerr << "Unknown exception!" << '\n'
                << "Aborting!" << '\n'
                << "----------------------------------------------------" << '\n';
      return 1;
    }

  return 0;
}
