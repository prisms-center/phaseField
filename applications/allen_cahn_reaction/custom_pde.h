// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>

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
    , m_mix(get_user_inputs().user_constants.get_double("m_mix"))
    , kappa(get_user_inputs().user_constants.get_double("kappa"))
    , c1(get_user_inputs().user_constants.get_double("c1"))
    , c2(get_user_inputs().user_constants.get_double("c2"))
    , l_int(get_user_inputs().user_constants.get_double("l_int"))
    , Iapp(get_user_inputs().user_constants.get_double("Iapp")) // applied current
    , Vapp(get_user_inputs().user_constants.get_double("Vapp")) // initial voltage guess
    , tol_Vapp(get_user_inputs().user_constants.get_double("tol_Vapp"))
    , max_newton(get_user_inputs().user_constants.get_int("max_newton"))
  {}

private:
  number m_mix;
  number kappa;
  number c1;
  number c2;
  number l_int;
  number Iapp;
  number Vapp;
  number tol_Vapp;
  int    max_newton;

  // reaction
  template <typename ValueType>
  [[nodiscard]] inline ValueType
  jbv(const ValueType &eta) const
  {
    return 2.0 * std::sinh(0.5 * eta);
  }

  template <typename ValueType>
  [[nodiscard]] inline ValueType
  djbv(const ValueType &eta) const
  {
    return std::cosh(0.5 * eta);
  }

  // derivative of the chemical free energy density
  [[nodiscard]] inline ScalarValue
  dfchem(const ScalarValue &c) const
  {
    return m_mix * (1.0 - 2.0 * c) + std::log(c / (1.0 - c));
  }

  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    if (index == 0) // c
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
        constexpr number rad[12] =
          {0.12, 0.14, 0.19, 0.16, 0.11, 0.12, 0.17, 0.15, 0.20, 0.10, 0.11, 0.14};
        number dist = 0.0;
        number sdf  = std::numeric_limits<number>::max();
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
        scalar_value = c1 + (c2 - c1) * (0.5 * (1.0 - std::tanh(sdf / (0.5 * l_int))));
      }
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    if (solve_block_id == 0) // explicit: c
      {
        ScalarValue c  = variable_list.template get_value<Scalar, OldOne>(0);
        ScalarValue mu = variable_list.template get_value<Scalar, OldOne>(1);

        ScalarValue R    = jbv(mu + ScalarValue(Vapp));
        ScalarValue eq_c = c + sim_timer.get_timestep() * R;

        variable_list.set_value_term(0, eq_c);
      }
    else if (solve_block_id == 1) // mu
      {
        ScalarValue c      = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad  c_grad = variable_list.template get_gradient<Scalar, Current>(0);

        variable_list.set_value_term(1, dfchem(c));
        variable_list.set_gradient_term(1, kappa * c_grad);
      }
    else if (solve_block_id == 2) // postprocess
      {
        ScalarValue mu = variable_list.template get_value<Scalar, Current>(1);

        variable_list.set_value_term(2, jbv(mu + ScalarValue(Vapp)));
      }
  }

  void
  post_solve_block([[maybe_unused]] SolveContext<dim, degree, number> &solve_context,
                   [[maybe_unused]] unsigned int solve_block_id) override
  {
    if (solve_block_id == 1) // mu
      {
        // impose galvanostatic condition
        const auto &dof_handler =
          solve_context.get_dof_manager().get_field_dof_handler(1);
        const auto &mu_solution =
          solve_context.get_solution_indexer().get_solution_vector(1);

        number V_guess = Vapp; // current guess
        number dV {};

        unsigned int i;
        for (i = 0; i < max_newton; i++)
          {
            number I1 = Integrator<dim, degree, number>::integrate_function(
              dof_handler,
              [this, V_guess](const number mu_q) -> number
              {
                return this->jbv(mu_q + V_guess);
              },
              mu_solution);
            number I2 = Integrator<dim, degree, number>::integrate_function(
              dof_handler,
              [this, V_guess](const number mu_q) -> number
              {
                return this->djbv(mu_q + V_guess);
              },
              mu_solution);

            dV = (Iapp - I1) / I2;
            V_guess += dV;

            if (std::abs(dV) < tol_Vapp)
              {
                break;
              }
          }
        Vapp = V_guess;

        Logger::instance() << LogFormatter::info("Applied voltage (Vapp): ") << Vapp
                           << ", |dV|: " << std::abs(dV) << ", niter:" << i + 1 << "\n";
      }
  }
};

PRISMS_PF_END_NAMESPACE
