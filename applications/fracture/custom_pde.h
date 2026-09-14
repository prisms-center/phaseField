// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/pde_operator_base.h>
#include <prismspf/utilities/mechanics.h>

#include <prismspf/config.h>

#include <cmath>
#include <numbers>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public PDEOperatorBase<dim, degree, number>
{
    public:
        using ScalarValue = dealii::VectorizedArray<number>;
        using ScalarGrad  = dealii::Tensor<1, dim, ScalarValue>;
        using VectorGrad  = dealii::Tensor<2, dim, ScalarValue>;
        using PDEOperatorBase<dim, degree, number>::get_user_inputs;
        using PDEOperatorBase<dim, degree, number>::get_pf_tools;

        CustomPDE(const UserInputParameters<dim> &_user_inputs, PhaseFieldTools<dim> &_pf_tools)
            : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
            , clength(get_user_inputs().user_constants.get_double("cracklength"))
            , Mn(get_user_inputs().user_constants.get_double("Mn"))
            , ell(get_user_inputs().user_constants.get_double("ell"))
            , Gc0(get_user_inputs().user_constants.get_double("Gc0"))
            , CIJ_base(get_user_inputs().user_constants.get_elasticity_tensor("CIJ_base"))
            , KI_nom(get_user_inputs().user_constants.get_double("KI_nom"))
            , vel_nom(get_user_inputs().user_constants.get_double("vel_nom"))
        {
            const auto &sd = get_user_inputs().spatial_discretization;
            dx = static_cast<number>(sd.rectangular_mesh.size[0])
                / static_cast<number>(sd.rectangular_mesh.subdivisions[0])
                / std::pow(number(2.0), static_cast<number>(sd.global_refinement));
            Ly = static_cast<number>(sd.rectangular_mesh.size[1]);
        }

    private:
        void
        set_initial_condition([[maybe_unused]] const unsigned int       &index,
                              [[maybe_unused]] const unsigned int       &component,
                              [[maybe_unused]] const dealii::Point<dim> &point,
                              [[maybe_unused]] number                   &scalar_value,
                              [[maybe_unused]] number &vector_component_value) const override
        {
            scalar_value           = 0.0;
            vector_component_value = 0.0;

            if (index == 0) // n: horizontal crack seed at mid-height, x < clength
            {
                if (std::abs(point[1] - (Ly / 2.0) + (0.5 * dx)) < dx && point[0] < clength) {
                    scalar_value = 1.0;
                }
            }
            if (index == 3 || index == 4) { // Ex, Gx = 1 everywhere (homogeneous material)
                scalar_value = 1.0;
            }
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
            if (index == 1)
            {
                constexpr double pi = std::numbers::pi;

                number x     = point[0] - (vel_nom * static_cast<number>(sim_timer.get_time())) - clength;
                number y     = point[1] - (Ly / 2.0) + (dx * 0.5);
                number r     = std::sqrt((x * x) + (y * y));
                number theta = std::atan2(y, x);
                number mu    = CIJ_base[dim][dim];
                number lam   = CIJ_base[0][0] - (2.0 * mu);
                number nu    = lam / (2.0 * (lam + mu));
                number kappa = 3.0 - 4.0 * nu;

                number val = 0.5 * (KI_nom / mu) * std::sqrt(0.5 * r / pi) * (kappa - std::cos(theta));
                if (component == 0) {
                    vector_component_value = val * std::cos(0.5 * theta);
                } else {
                    vector_component_value = val * std::sin(0.5 * theta);
                }
            }
        }

        void
        compute_rhs(FieldContainer<dim, degree, number> &variable_list,
                    const SimulationTimer               &sim_timer,
                    unsigned int                         solve_block_id) const override
        {
            if (solve_block_id == 0) // explicit n update
            {
                ScalarValue n    = variable_list.template get_value<Scalar, OldOne>(0);
                ScalarValue dndt = variable_list.template get_value<Scalar, OldOne>(2);
                const number dt  = static_cast<number>(sim_timer.get_timestep());

                for (unsigned int j = 0; j < dndt.size(); ++j)
                {
                    if (dndt[j] > 0.0) {
                        dndt[j] = 0.0;
                    }
                    if (n[j] - (dndt[j] * dt) > 1.0) {
                        dndt[j] = (n[j] - 1.0) / dt;
                    }
                }
                variable_list.set_value_term(0, n - (dt * dndt));
            }

            else if (solve_block_id == 3) // auxiliary dndt
            {
                ScalarValue n  = variable_list.template get_value<Scalar, Current>(0);
                ScalarGrad  nx = variable_list.template get_gradient<Scalar, Current>(0);
                VectorGrad  ux = variable_list.template get_symmetric_gradient<Vector, Current>(1);
                ScalarValue Ex = variable_list.template get_value<Scalar, Current>(3);
                ScalarValue Gx = variable_list.template get_value<Scalar, Current>(4);

                dealii::Tensor<2, Mechanics::voigt_tensor_size<dim>, ScalarValue> C = CIJ_base * Ex;
                VectorGrad stress;
                Mechanics::compute_stress<dim, ScalarValue>(C, ux, stress);
                ScalarValue psi = 0.0;
                for (unsigned int i = 0; i < dim; ++i) {
                    for (unsigned int j = 0; j < dim; ++j) {
                        psi += 0.5 * stress[i][j] * ux[i][j];
                    }
                }

                variable_list.set_value_term(
                    2, (2.0 * (n - 1.0) * psi + Gc0 * Gx * 3.0 / 8.0 / ell) * Mn);
                variable_list.set_gradient_term(2, ell * nx * Gc0 * Gx * (3.0 / 4.0) * Mn);
            }

            else if (solve_block_id == 2) // no body force; analytic BC applied via set_dirichlet
            {
                variable_list.set_gradient_term(1, VectorGrad());
            }

            else if (solve_block_id == 4) // postprocessing
            {
                ScalarValue n  = variable_list.template get_value<Scalar, Current>(0);
                ScalarGrad  nx = variable_list.template get_gradient<Scalar, Current>(0);
                VectorGrad  ux = variable_list.template get_symmetric_gradient<Vector, Current>(1);
                ScalarValue Ex = variable_list.template get_value<Scalar, Current>(3);
                ScalarValue Gx = variable_list.template get_value<Scalar, Current>(4);

                ScalarValue f_int =
                    (Gc0 * n * Gx * (3.0 / 8.0) / ell) + (Gc0 * Gx * 3.0 / 8.0 * ell * nx * nx);

                dealii::Tensor<2, Mechanics::voigt_tensor_size<dim>, ScalarValue> C_deg =
                    CIJ_base * Ex * (1.0 - 2.0 * n + n * n);
                VectorGrad stress;
                Mechanics::compute_stress<dim, ScalarValue>(C_deg, ux, stress);
                ScalarValue f_el = 0.0;
                for (unsigned int i = 0; i < dim; ++i) {
                    for (unsigned int j = 0; j < dim; ++j) {
                        f_el += 0.5 * stress[i][j] * ux[i][j];
                    }
                }

                variable_list.set_value_term(5, f_el + f_int); // f_tot
                variable_list.set_value_term(6, stress[0][0]); // s11
                variable_list.set_value_term(7, stress[0][1]); // s12
                variable_list.set_value_term(8, stress[1][1]); // s22
                variable_list.set_value_term(9, f_int);
                variable_list.set_value_term(10, f_el);
            }
        }

        void
        compute_lhs(FieldContainer<dim, degree, number>    &variable_list,
                    [[maybe_unused]] const SimulationTimer &sim_timer,
                    unsigned int                            solve_block_id) const override
        {
            if (solve_block_id == 2) // degraded stiffness tangent for u
            {
                ScalarValue n        = variable_list.template get_value<Scalar, Current>(0);
                VectorGrad  ux_trial = variable_list.template get_symmetric_gradient<Vector, LHS>(1);
                ScalarValue Ex       = variable_list.template get_value<Scalar, Current>(3);

                dealii::Tensor<2, Mechanics::voigt_tensor_size<dim>, ScalarValue> C_deg =
                    CIJ_base * Ex * (1.0 - 2.0 * n + n * n);
                VectorGrad stress;
                Mechanics::compute_stress<dim, ScalarValue>(C_deg, ux_trial, stress);
                variable_list.set_gradient_term(1, stress);
            }
        }

        // ---- member variables ----

        number clength;
        number Mn;
        number ell;
        number Gc0;
        dealii::Tensor<2, Mechanics::voigt_tensor_size<dim>, number> CIJ_base;
        number KI_nom;
        number vel_nom;
        number dx = 0.0;
        number Ly = 0.0;

};

PRISMS_PF_END_NAMESPACE
