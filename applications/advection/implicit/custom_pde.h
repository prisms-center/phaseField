// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/function_signed_distance.h>

#include <prismspf/core/pde_operator_base.h>

#include <prismspf/utilities/utilities.h>

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

  CustomPDE(const UserInputParameters<dim> &_user_inputs, PhaseFieldTools<dim> &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , W(_user_inputs.user_constants.get_double("W"))
  {
    // Compute the center of the domain, promoting the doubles to vectorized arrays
    // NOTE: This assumes a rectangular mesh
    for (unsigned int d = 0; d < dim; ++d)
      {
        domain_center[d] =
          _user_inputs.spatial_discretization.rectangular_mesh.size[d] / 2.0;
      }
  };

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    using std::tanh;

    // Zalesak's disk level-set
    // TODO: Refactor these
    const double             radius       = 15.0;
    const double             notch_width  = 5.0;
    const double             notch_height = 25.0;
    const dealii::Point<dim> center(50.0, 75.0);

    dealii::Functions::SignedDistance::ZalesakDisk<dim> zalesak_disk(center,
                                                                     radius,
                                                                     notch_width,
                                                                     notch_height);
    // The distance is just the level-set value
    const double distance = zalesak_disk.value(point);

    // Apply tanh
    scalar_value = 0.5 * (1.0 - tanh(distance / (std::numbers::sqrt2 * W)));
  };

  void
  set_dirichlet([[maybe_unused]] const unsigned int       &index,
                [[maybe_unused]] const unsigned int       &boundary_id,
                [[maybe_unused]] const unsigned int       &component,
                [[maybe_unused]] const dealii::Point<dim> &point,
                [[maybe_unused]] number                   &scalar_value,
                [[maybe_unused]] number &vector_component_value) const override {};

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_group_id) const override
  {
    ScalarValue dt = sim_timer.get_timestep();

    if (solve_group_id == 0)
      {
        ScalarValue u = variable_list.template get_value<Scalar, OldOne>(0);
        ScalarValue h = variable_list.get_element_volume();
        dealii::Point<dim, ScalarValue> p = variable_list.get_q_point_location();

        // Compute the rotational velocity based on the distance from the center of the
        // domain
        // TODO: This assumes 2D
        ScalarGrad v;
        v[0] = -(p[1] - domain_center[1]);
        v[1] = p[0] - domain_center[0];

        // Compute the stabilization parameter
        ScalarValue tau = stabilization_parameter<dim, degree>(dt, h, v);

        // Compute the RHS of residual for SUPG stabilization
        ScalarValue R = -u;

        variable_list.set_value_term(0, u);
        variable_list.set_gradient_term(0, tau * R * v);
      }
  };

  void
  compute_lhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solver_id) const override
  {
    ScalarValue dt = sim_timer.get_timestep();

    if (solver_id == 0)
      {
        ScalarValue u      = variable_list.template get_value<Scalar, LHS>(0);
        ScalarGrad  u_grad = variable_list.template get_gradient<Scalar, LHS>(0);
        ScalarValue h      = variable_list.get_element_volume();
        dealii::Point<dim, ScalarValue> p = variable_list.get_q_point_location();

        // Compute the rotational velocity based on the distance from the center of the
        // domain
        // TODO: This assumes 2D
        ScalarGrad v;
        v[0] = -(p[1] - domain_center[1]);
        v[1] = p[0] - domain_center[0];

        // Compute the stabilization parameter
        ScalarValue tau = stabilization_parameter<dim, degree>(dt, h, v);

        // Compute the LHS of residual for SUPG stabilization
        ScalarValue R = u - dt * v * u_grad;

        variable_list.set_value_term(0, u);
        variable_list.set_gradient_term(0, -v * dt * u - tau * R * v);
      }
  };

  number                          W;
  dealii::Point<dim, ScalarValue> domain_center;
};

PRISMS_PF_END_NAMESPACE
