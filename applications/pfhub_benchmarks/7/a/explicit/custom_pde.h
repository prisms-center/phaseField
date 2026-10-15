// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>

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
  CustomPDE(const UserInputParameters<dim> &_user_inputs, PhaseFieldTools<dim> &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , kappa(0.0004)
    , A_1(0.0075)
    , B_1(8.0*std::numbers::pi)
    , A_2(0.03)
    , B_2(22.0*std::numbers::pi)
    , C_2(0.0625*std::numbers::pi)
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    scalar_value = eta0(point[0], point[1]);
  }

void
  set_dirichlet([[maybe_unused]] const unsigned int       &index,
                [[maybe_unused]] const unsigned int       &boundary_id,
                [[maybe_unused]] const unsigned int       &component,
                [[maybe_unused]] const dealii::Point<dim> &point,
                [[maybe_unused]] const SimulationTimer    &sim_timer,
                [[maybe_unused]] number                   &scalar_value,
                [[maybe_unused]] number &vector_component_value) const override
  {
    if (boundary_id == RectangularMesh<dim>::Boundary::Bottom)
      {
        scalar_value = 1.0;
      }
    if (boundary_id == RectangularMesh<dim>::Boundary::Top)
      {
        scalar_value = 0.0;
      }
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_block_id) const override
  {
    dealii::Point<dim, ScalarValue> point = variable_list.get_q_point_location();
    if (solve_block_id == 1) // explicit n
      {
        ScalarValue n_val  = variable_list.template get_value<Scalar, OldOne>(0);
        ScalarGrad  n_grad = variable_list.template get_gradient<Scalar, OldOne>(0);
        ScalarValue f_well = 4.0 * n_val * (n_val - 1.0) * (n_val - 0.5);
        ScalarValue S_val = S(point[0], point[1], sim_timer.get_time());

        ScalarValue eq_n   = n_val - sim_timer.get_timestep() * (f_well - S_val);
        ScalarGrad  eqx_n  = -sim_timer.get_timestep() * kappa * n_grad;

        variable_list.set_value_term(0, eq_n);
        variable_list.set_gradient_term(0, eqx_n);
      }
    else if (solve_block_id == 2) // postprocess
      {
        ScalarValue n_val  = variable_list.template get_value<Scalar, Current>(0);
        ScalarValue n_sol = eta(point[0], point[1], sim_timer.get_time());
        ScalarValue error = n_val - n_sol;

        variable_list.set_value_term(1, n_sol);
        variable_list.set_value_term(2, error);
      }
  }

  number kappa;
  number A_1;
  number B_1;
  number A_2;
  number B_2;
  number C_2;

  // Custom functions for PFHub Benchmark 7
  // Generated with SymPy using the included script: generate_cpp_equations.py
  [[nodiscard]] double eta0(double x, double y) const
  {
    return (1.0/2.0)*std::tanh((1.0/8.0)*std::numbers::sqrt2*(4*A_2*std::sin(B_2*x) - 4.0*y + 1.0)/std::sqrt(kappa)) + 0.5;
  }

  ScalarValue eta(ScalarValue x, ScalarValue y, double t) const
  {
    return (1.0/2.0)*std::tanh((1.0/8.0)*std::numbers::sqrt2*(4*A_1*t*std::sin(B_1*x) + 4.0*A_2*std::sin(B_2*x + C_2*t) - 4.0*y + 1.0)/std::sqrt(kappa)) + 0.5;
  }

  ScalarValue alpha(ScalarValue x, double t) const
  {
    return A_1*t*std::sin(B_1*x) + A_2*std::sin(B_2*x + C_2*t) + 0.25;
  }

  ScalarValue sech(ScalarValue arg) const
  {
    return 1.0 / std::cosh(arg);
  }

  ScalarValue S(ScalarValue x, ScalarValue y, double t) const
  {
    return (-2.0*std::sqrt(kappa)*std::pow((A_1*B_1*t*std::cos(B_1*x) + A_2*B_2*std::cos(B_2*x + C_2*t)), 2.0)*std::tanh(std::numbers::sqrt2*(-A_1*t*std::sin(B_1*x) - A_2*std::sin(B_2*x + C_2*t) + y - 0.25)/(2.0*std::sqrt(kappa))) - std::numbers::sqrt2*kappa*(-A_1*B_1*B_1*t*std::sin(B_1*x) - A_2*B_2*B_2*std::sin(B_2*x + C_2*t)) + std::numbers::sqrt2*(A_1*std::sin(B_1*x) + A_2*C_2*std::cos(B_2*x + C_2*t)))*std::pow(sech(std::numbers::sqrt2*(-A_1*t*std::sin(B_1*x) - A_2*std::sin(B_2*x + C_2*t) + y - 0.25)/(2.0*std::sqrt(kappa))), 2.0)/(4.0*std::sqrt(kappa));
  }
};

PRISMS_PF_END_NAMESPACE
