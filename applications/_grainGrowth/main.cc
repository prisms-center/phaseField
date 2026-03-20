// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>

#include <prismspf/core/dependencies.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/grid_refiner_criterion.h>
#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/phase_field_tools.h>
#include <prismspf/core/problem.h>
#include <prismspf/core/simulation_timer.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/input_file_reader.h>
#include <prismspf/user_inputs/miscellaneous_parameters.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <random>

using namespace prisms;

constexpr unsigned int num_ops = 5;

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
    , m(get_user_inputs().user_constants.get_model_constant_double("m"))
    , alpha(get_user_inputs().user_constants.get_model_constant_double("alpha"))
    , kappa(get_user_inputs().user_constants.get_model_constant_double("kappa"))
    , L(get_user_inputs().user_constants.get_model_constant_double("L"))
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    using dealii::Point;
    const dealii::Tensor<1, dim> &mesh_size =
      get_user_inputs().spatial_discretization.rectangular_mesh.size;
    if (index < 5)
      {
        std::vector<Point<dim>> center_list;

        // The big grains
        {
          Point<dim> center(0.2, 0.15);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.25, 0.7);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.5, 0.5);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.6, 0.85);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.85, 0.35);
          center_list.push_back(center);
        }

        // The medium grains
        {
          Point<dim> center(0.08, 0.92);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.75, 0.6);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.75, 0.1);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.2, 0.45);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.85, 0.85);
          center_list.push_back(center);
        }

        // The small grains
        {
          Point<dim> center(0.55, 0.05);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.1, 0.35);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.95, 0.65);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.9, 0.15);
          center_list.push_back(center);
        }
        {
          Point<dim> center(0.45, 0.25);
          center_list.push_back(center);
        }

        std::vector<double> rad = {0.14,
                                   0.14,
                                   0.14,
                                   0.14,
                                   0.14,
                                   0.08,
                                   0.08,
                                   0.08,
                                   0.08,
                                   0.08,
                                   0.05,
                                   0.05,
                                   0.05,
                                   0.05,
                                   0.05};

        double dist  = 0.0;
        scalar_value = 0;

        for (unsigned int dir = 0; dir < dim; dir++)
          {
            dist += (point[dir] - center_list[index][dir] * mesh_size[dir]) *
                    (point[dir] - center_list[index][dir] * mesh_size[dir]);
          }
        dist = std::sqrt(dist);

        scalar_value += 0.5 * (1.0 - std::tanh((dist - rad[index] * mesh_size[0]) / 0.5));

        dist = 0.0;
        for (unsigned int dir = 0; dir < dim; dir++)
          {
            dist += (point[dir] - center_list[index + 5][dir] * mesh_size[dir]) *
                    (point[dir] - center_list[index + 5][dir] * mesh_size[dir]);
          }
        dist = std::sqrt(dist);

        scalar_value +=
          0.5 * (1.0 - std::tanh((dist - rad[index + 5] * mesh_size[0]) / 0.5));

        dist = 0.0;
        for (unsigned int dir = 0; dir < dim; dir++)
          {
            dist += (point[dir] - center_list[index + 10][dir] * mesh_size[dir]) *
                    (point[dir] - center_list[index + 10][dir] * mesh_size[dir]);
          }
        dist = std::sqrt(dist);

        scalar_value +=
          0.5 * (1.0 - std::tanh((dist - rad[index + 10] * mesh_size[0]) / 0.5));
      }
    else
      {
        scalar_value = 0.0;
      }
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_group_id) const override
  {
    double dt = sim_timer.get_timestep();
    if (solve_group_id == 0) // grains
      {
        dealii::AlignedVector<ScalarValue> n_val(num_ops);
        dealii::AlignedVector<ScalarGrad>  n_grad(num_ops);
        dealii::AlignedVector<ScalarValue> dfdn_val(num_ops);
        dealii::AlignedVector<ScalarGrad>  dfdn_vec(num_ops);
        ScalarValue                        sum_n_sq;
        for (unsigned int op_id = 0; op_id < num_ops; op_id++)
          {
            const ScalarValue n =
              variable_list.template get_value<TensorRank::Scalar, OldOne>(op_id);
            const ScalarValue nsq = n * n;
            n_val[op_id]          = n;
            sum_n_sq += nsq;
            dfdn_val[op_id] = nsq * n - n;
            dfdn_vec[op_id] =
              -kappa *
              variable_list.template get_gradient<TensorRank::Scalar, OldOne>(op_id);
          }
        for (unsigned int op_id = 0; op_id < num_ops; op_id++)
          {
            const ScalarValue &n = n_val[op_id];
            dfdn_val[op_id] += 2.0 * alpha * n * (sum_n_sq - n * n);
            dfdn_val[op_id] *= m;
          }
        for (unsigned int op_id = 0; op_id < num_ops; op_id++)
          {
            variable_list.set_value_term(op_id, n_val[op_id] + dt * -L * dfdn_val[op_id]);
            variable_list.set_gradient_term(op_id, -dt * -L * dfdn_vec[op_id]);
          }
      }
    else if (solve_group_id == 1) // pp
      {
        ScalarValue sum_n_sq;
        for (unsigned int op_id = 0; op_id < num_ops; op_id++)
          {
            const ScalarValue n =
              variable_list.template get_value<TensorRank::Scalar, Normal>(op_id);
            sum_n_sq += n * n;
          }
        variable_list.set_value_term(num_ops, 1.0 * (1.0 - sum_n_sq));
      }
  }

  number m;
  number alpha;
  number kappa;
  number L;
};

int
main(int argc, char *argv[])
{
  // Initialize MPI
  dealii::Utilities::MPI::MPI_InitFinalize
    mpi_init(argc, argv, dealii::numbers::invalid_unsigned_int);

  // Restrict deal.II console printing
  dealii::deallog.depth_console(0);

  // Parse the command line options (if any) to get the name of the input file
  ParseCMDOptions cli_options(argc, argv);

  constexpr unsigned int dim    = 2;
  constexpr unsigned int degree = 2;

  // Declare grains
  std::vector<FieldAttributes> fields(num_ops + 1);
  for (unsigned int op_id = 0; op_id < num_ops; op_id++)
    {
      fields[op_id] = FieldAttributes("n" + std::to_string(op_id));
    }
  fields[num_ops] = FieldAttributes("gb");

  // Set up grain solver
  SolveGroup grain_solve;
  grain_solve.id           = 0;
  grain_solve.solve_type   = Explicit;
  grain_solve.solve_timing = Initialized;
  const Dependency grain_dependency(EvalFlags::nothing,
                                    EvalFlags::nothing,
                                    {EvalFlags::values | EvalFlags::gradients});
  for (unsigned int op_id = 0; op_id < num_ops; op_id++)
    {
      grain_solve.field_indices.insert(op_id);
      grain_solve.dependencies_rhs[op_id] = grain_dependency;
    }

  // Set up postprocess solver
  SolveGroup pp_group;
  pp_group.id            = 1;
  pp_group.solve_type    = Explicit;
  pp_group.solve_timing  = PostProcess;
  pp_group.field_indices = {num_ops};
  const Dependency pp_dependency(EvalFlags::values);
  for (unsigned int op_id = 0; op_id < num_ops; op_id++)
    {
      pp_group.dependencies_rhs[op_id] = pp_dependency;
    }

  const std::vector<SolveGroup> solve_groups({grain_solve, pp_group});

  UserInputParameters<dim> user_inputs(cli_options.get_parameters_filename());
  // set refinement criterion for all grains
  for (unsigned int op_id = 0; op_id < num_ops; op_id++)
    {
      RefinementCriterion refine_criterion(RefinementFlags::Value, 0.001, 0.999);
      user_inputs.spatial_discretization.refinement_criteria[fields[op_id].name] =
        refine_criterion;
    }
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
