// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>

#include <prismspf/solvers/solver_base.h>

#include <prismspf/utilities/integrator.h>

#include <cmath>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public PDEOperatorBase<dim, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, ScalarValue>;
  using ScalarHess  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;
  using VectorGrad  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorHess  = dealii::Tensor<3, dim, ScalarValue>;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  /**
   * @brief Constructor.
   */
  explicit CustomPDE(const UserInputParameters<dim> &_user_inputs,
                     PhaseFieldTools<dim>           &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , mobility(get_user_inputs().user_constants.get_double("mobility"))
    , kappa(get_user_inputs().user_constants.get_double("kappa"))
  {}

  void
  set_lambda(number v)
  {
    lambda = v;
  }

private:
  number mobility;
  number kappa;
  number lambda;

  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    const dealii::Tensor<1, dim> &mesh_size =
      get_user_inputs().spatial_discretization.rectangular_mesh.size;
    constexpr number center[12][3] = {
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
    constexpr number rad[12] = {12, 14, 19, 16, 11, 12, 17, 15, 20, 10, 11, 14};

    number dist = 0.0;
    number sdf  = std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < 12; i++)
      {
        dist = 0.0;
        for (unsigned int dir = 0; dir < dim; dir++)
          {
            const number comp_diff = point[dir] - center[i][dir] * mesh_size[dir];
            dist += comp_diff * comp_diff;
          }
        dist = std::sqrt(dist) - rad[i];
        sdf  = std::min(sdf, dist);
      }
    scalar_value = 0.5 * (1.0 - std::tanh(sdf / 1.5));
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_block_id) const override
  {
    if (solve_block_id == 0) // explicit c
      {
        ScalarValue c  = variable_list.template get_value<Scalar, OldOne>(0);
        ScalarValue mu = variable_list.template get_value<Scalar, OldOne>(1);

        ScalarValue eq_c = c - sim_timer.get_timestep() * mobility * (mu - lambda);

        variable_list.set_value_term(0, eq_c);
      }
    else if (solve_block_id == 1) // mu
      {
        ScalarValue c       = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad  c_grad  = variable_list.template get_gradient<Scalar, Current>(0);
        ScalarValue df_well = 4.0 * c * (c - 1.0) * (c - 0.5);

        variable_list.set_value_term(1, df_well);
        variable_list.set_gradient_term(1, kappa * c_grad);
      }
    else if (solve_block_id == 2) // custom
      {
      }
    else if (solve_block_id == 3) // postprocess
      {
        ScalarValue c      = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad  c_grad = variable_list.template get_gradient<Scalar, Current>(0);
        ScalarValue f_chem = c * c * (1.0 - c) * (1.0 - c);
        ScalarValue f_grad = 0.5 * kappa * c_grad.norm_square();

        variable_list.set_value_term(3, c_grad.norm());
        variable_list.set_value_term(4, f_chem + f_grad);
      }
  }
};

template <unsigned int dim, unsigned int degree, typename number>
class MyCustomSolver : public SolverBase<dim, degree, number>
{
private:
  CustomPDE<dim, degree, number> &custom_pde;

public:
  MyCustomSolver(const SolveBlock                        &solve_block,
                 const SolveContext<dim, degree, number> &solve_context,
                 CustomPDE<dim, degree, number>          &_custom_pde)
    : SolverBase<dim, degree, number>(solve_block, solve_context)
    , custom_pde(_custom_pde)
  {}

  void
  init(const std::list<DependencyMap> &all_dependency_sets) override
  {
    SolverBase<dim, degree, number>::init(all_dependency_sets);
  }

  void
  solve_level(unsigned int relative_level) override
  {
    number mu_int = Integrator<dim, degree, number>::template integrate<0>(
      this->solve_context->get_dof_manager().get_field_dof_handler(1),
      this->solve_context->get_solution_indexer().get_solution_vector(1))[0];

    const dealii::Tensor<1, dim> &mesh_size =
      this->solve_context->get_user_inputs().spatial_discretization.rectangular_mesh.size;
    number V = mesh_size[0] * mesh_size[1];

    custom_pde.set_lambda(mu_int / V);
  }

  void
  update() override
  {}
};

PRISMS_PF_END_NAMESPACE
