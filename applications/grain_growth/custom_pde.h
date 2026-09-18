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

  CustomPDE(const UserInputParameters<dim> &_user_inputs, PhaseFieldTools<dim> &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , m_well(get_user_inputs().user_constants.get_double("m_well"))
    , kappa(get_user_inputs().user_constants.get_double("kappa"))
    , alpha(get_user_inputs().user_constants.get_double("alpha"))
    , kinetic_coef(get_user_inputs().user_constants.get_double("kinetic_coef"))
  {};

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    // For this example, the initial grain structure is given by a set of
    // overlapping circles defined by a hyperbolic tangent function. The
    // center of each circle is given by "center_list" and the radius by
    // "radius"

    if (index < 5)
      {
        const dealii::Tensor<1, dim> &mesh_size =
          get_user_inputs().spatial_discretization.rectangular_mesh.size;

        std::vector<dealii::Point<dim>> center_list;
        std::vector<number> radius = {
          0.14, 0.14, 0.14, 0.14, 0.14,
          0.08, 0.08, 0.08, 0.08, 0.08,
          0.05, 0.05, 0.05, 0.05, 0.05
        };

        // Big grains
        center_list.push_back(dealii::Point<dim>(0.2, 0.15));
        center_list.push_back(dealii::Point<dim>(0.25, 0.7));
        center_list.push_back(dealii::Point<dim>(0.5, 0.5));
        center_list.push_back(dealii::Point<dim>(0.6, 0.85));
        center_list.push_back(dealii::Point<dim>(0.85, 0.35));

        // Medium grains
        center_list.push_back(dealii::Point<dim>(0.08, 0.92));
        center_list.push_back(dealii::Point<dim>(0.75, 0.6));
        center_list.push_back(dealii::Point<dim>(0.75, 0.1));
        center_list.push_back(dealii::Point<dim>(0.2, 0.45));
        center_list.push_back(dealii::Point<dim>(0.85, 0.85));

        // Small grains
        center_list.push_back(dealii::Point<dim>(0.55, 0.05));
        center_list.push_back(dealii::Point<dim>(0.1, 0.35));
        center_list.push_back(dealii::Point<dim>(0.95, 0.65));
        center_list.push_back(dealii::Point<dim>(0.9, 0.15));
        center_list.push_back(dealii::Point<dim>(0.45, 0.25));

        scalar_value = 0.0;

        // The big grains
        double dist = 0.0;
        for (unsigned int dir = 0; dir < dim; dir++)
          {
            dist += (point[dir] - center_list[index][dir] * mesh_size[dir])
              * (point[dir] - center_list[index][dir] * mesh_size[dir]);
          }
        dist = std::sqrt(dist);
        scalar_value += 0.5*(1.0 - std::tanh((dist - radius[index] * mesh_size[0]) / 0.5));

        // The medium grains
        dist = 0.0;
        for (unsigned int dir = 0; dir < dim; dir++)
          {
            dist += (point[dir] - center_list[index + 5][dir] * mesh_size[dir])
              * (point[dir] - center_list[index + 5][dir] * mesh_size[dir]);
          }
        dist = std::sqrt(dist);
        scalar_value += 0.5*(1.0 - std::tanh((dist - radius[index + 5] * mesh_size[0]) / 0.5));

        // The small grains
        dist = 0.0;
        for (unsigned int dir = 0; dir < dim; dir++)
          {
            dist += (point[dir] - center_list[index + 10][dir] * mesh_size[dir])
              * (point[dir] - center_list[index + 10][dir] * mesh_size[dir]);
          }
        dist = std::sqrt(dist);
        scalar_value += 0.5*(1.0 - std::tanh((dist - radius[index + 10] * mesh_size[0]) / 0.5));
      }
    else
      {
        scalar_value = 0.0;
      }
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override 
  {
    if (solve_block_id == 1) // Explicit solve for the order parameters
      {
        ScalarValue f_multiwell;
        ScalarValue ni;
        ScalarValue nj;
        ScalarGrad nix;

        std::vector<ScalarValue> value_terms(number_of_fields);
        std::vector<ScalarGrad> gradient_terms(number_of_fields);

        // Calculate the evolution equations
        for (unsigned int i = 0; i < number_of_fields; i++)
          {
            ni  = variable_list.template get_value<Scalar, OldOne>(i);
            nix = variable_list.template get_gradient<Scalar, OldOne>(i);
            f_multiwell = -ni + ni * ni * ni;
            for (unsigned int j = 0; j < number_of_fields; j++)
              {
                if (i != j)
                  {
                    nj = variable_list.template get_value<Scalar, OldOne>(j);
                    f_multiwell += 2.0 * alpha * ni * nj * nj;
                  }
              }
            value_terms[i] = ni - sim_timer.get_timestep() * kinetic_coef * m_well * f_multiwell;
            gradient_terms[i] = -sim_timer.get_timestep() * kinetic_coef * kappa * nix;
          }

        // Submit the terms
        for (unsigned int i = 0; i < number_of_fields; i++)
          {
            variable_list.set_value_term(i, value_terms[i]);
            variable_list.set_gradient_term(i, gradient_terms[i]);
          }
      }
    else if (solve_block_id == 2) // Postprocessing
      {
        ScalarValue f_total = 0.0;
        ScalarValue sum2op = 0.0;
        ScalarValue ni;
        ScalarValue nj;
        ScalarGrad nix;

        ScalarValue max_op_id;
        ScalarValue max_op_value;

        for (unsigned int i = 0; i < number_of_fields; i++)
          {
            ni  = variable_list.template get_value<Scalar, Current>(i);
            nix = variable_list.template get_gradient<Scalar, Current>(i);
            f_total += m_well * (0.25 + (ni * ni * ni * ni / 4.0) - (ni * ni / 2.0));
            for (unsigned int j = 0; j < number_of_fields; j++)
              {
                if (i != j)
                  {
                    nj = variable_list.template get_value<Scalar, Current>(j);
                    f_total += m_well * alpha * ni * ni * nj * nj;
                    sum2op += ni * nj;
                  }
              }
            f_total += (kappa / 2.0) * nix.norm_square();

            for (unsigned int v = 0; v < ScalarValue::size(); v++)
              {
                if (ni[v] > max_op_value[v])
                  {
                    max_op_value[v] = ni[v];
                    max_op_id[v] = i;
                  }
              }
          }

        variable_list.set_value_term(number_of_fields, sum2op);
        variable_list.set_value_term(number_of_fields + 1, f_total);
        variable_list.set_value_term(number_of_fields + 2, max_op_id);
      }
  }

  number m_well;
  number kappa;
  number alpha;
  number kinetic_coef;
  
  static constexpr unsigned int number_of_fields = 6;
};

PRISMS_PF_END_NAMESPACE
