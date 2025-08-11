// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>

#include "catch.hpp"

#include <string>

/**
 * This unit test looks at variable_attributes.h and variable_attribute_loader.h and how
 * the parser determines the field sovle type of each field. In turn, the determines the
 * solver and other information such as nonlinearity.
 */
TEST_CASE("Field solve types")
{
  /**
   * Two nonlinear equations based on the steady-state version of the Cahn-Hilliard
   * equation Δ(u^3-u-γΔu)=0. In the weak form with γ=1 this becomes ∇w∇(u^3-u-Δu)=0.
   * This is a third order equation, so we solve the inner part with a coupled auxiliary
   * variable. ∇w∇x=0 and x = w(u^3-u)+(∇w∇u).
   *
   * Note that this scenario requires unintuitive selection of nonlinear criterion due to
   * the separation of variables.
   */
  SECTION("One time-independent and one auxiliary")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public variableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      loadVariableAttributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, SCALAR);
        set_variable_equation_type(0, TIME_INDEPENDENT);

        set_dependencies_value_term_LHS(0, "");
        set_dependencies_gradient_term_LHS(0, "");
        set_dependencies_value_term_RHS(0, "grad(eta)");
        set_dependencies_gradient_term_RHS(0, "");

        set_variable_name(1, "eta");
        set_variable_type(1, SCALAR);
        set_variable_equation_type(1, AUXILIARY);

        set_dependencies_value_term_LHS(1, "");
        set_dependencies_gradient_term_LHS(1, "");
        set_dependencies_value_term_RHS(1, "phi");
        set_dependencies_gradient_term_RHS(1, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    attributes.init_variable_attributes();
    AttributesList variables = attributes.get_var_attributes();

    REQUIRE(variables.size() == 2);
    for (unsigned int index : {0, 1})
      {
        REQUIRE(variables.at(index).field_solve_type ==
                fieldSolveType::NONEXPLICIT_CO_NONLINEAR);
      }
  }

  /**
   * Two nonlinear equations and two linear equations that are all independent of one
   * another. Note that the switching of vector and scalar field is intentional.
   *
   * The first equation is a steady-state version of the Burgers' equation  u∇u = vΔu.
   * In the weak form with v=1 this becomes w(u+δu)∇u + w(∇u+δ∇u)u = ∇w(∇u+δ∇u).
   *
   * The second equation is the minimal surface equation -∇⋅(∇u/√(1+|∇u|^2))=0. In the
   * weak form this becomes
   * ∇w⋅(∇u/√(1+|∇u|^2))+∇w⋅(δ∇u/√(1+|∇u|^2))-∇w⋅((∇u⋅δ∇u)∇u/√(1+|∇u|^2))=0.
   *
   * The third equation is a simple Poisson equation Δu=0.
   *
   * The fourth equation is the steady-state version of the Allen-Cahn equation
   * u^3-u-γΔu=0. In the weak form with γ=1 this becomes w(u^3+δu^3)+w(u-δu)+∇w(∇u+δ∇u)=0.
   */
  SECTION("Four time-independent")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public variableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      loadVariableAttributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, VECTOR);
        set_variable_equation_type(0, TIME_INDEPENDENT);

        set_dependencies_value_term_LHS(0,
                                        "change(phi), grad(change(phi)), phi, grad(phi)");
        set_dependencies_gradient_term_LHS(0, "grad(change(phi))");
        set_dependencies_value_term_RHS(0, "phi, grad(phi)");
        set_dependencies_gradient_term_RHS(0, "grad(phi)");

        set_variable_name(1, "eta");
        set_variable_type(1, SCALAR);
        set_variable_equation_type(1, TIME_INDEPENDENT);

        set_dependencies_value_term_LHS(1, "");
        set_dependencies_gradient_term_LHS(1, "grad(change(eta)), grad(eta)");
        set_dependencies_value_term_RHS(1, "");
        set_dependencies_gradient_term_RHS(1, "grad(eta)");

        set_variable_name(2, "beta");
        set_variable_type(2, VECTOR);
        set_variable_equation_type(2, TIME_INDEPENDENT);

        set_dependencies_value_term_LHS(2, "");
        set_dependencies_gradient_term_LHS(2, "grad(change(beta))");
        set_dependencies_value_term_RHS(2, "");
        set_dependencies_gradient_term_RHS(2, "grad(beta)");

        set_variable_name(3, "alpha");
        set_variable_type(3, SCALAR);
        set_variable_equation_type(3, TIME_INDEPENDENT);

        set_dependencies_value_term_LHS(3, "change(alpha)");
        set_dependencies_gradient_term_LHS(3, "grad(change(alpha))");
        set_dependencies_value_term_RHS(3, "alpha");
        set_dependencies_gradient_term_RHS(3, "grad(alpha)");
      }
    };

    testVariableAttributeLoader attributes;
    attributes.init_variable_attributes();
    AttributesList variables = attributes.get_var_attributes();

    REQUIRE(variables.size() == 4);
    for (unsigned int index : {0, 1})
      {
        REQUIRE(variables.at(index).field_solve_type ==
                fieldSolveType::NONEXPLICIT_SELF_NONLINEAR);
      }
    for (unsigned int index : {2, 3})
      {
        REQUIRE(variables.at(index).field_solve_type ==
                fieldSolveType::NONEXPLICIT_LINEAR);
      }
  }

  /**
   * Four explicit equations that are independent of one another.
   */
  SECTION("Four explicit")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public variableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      loadVariableAttributes() override
      {
        for (unsigned int index : {0, 1, 2, 3})
          {
            set_variable_name(index, "n" + std::to_string(index));
            set_variable_type(index, SCALAR);
            set_variable_equation_type(index, EXPLICIT_TIME_DEPENDENT);
            set_dependencies_value_term_RHS(index, "n" + std::to_string(index));
            set_dependencies_gradient_term_RHS(index,
                                               "grad(n" + std::to_string(index) + ")");
          }
      }
    };

    testVariableAttributeLoader attributes;
    attributes.init_variable_attributes();
    AttributesList variables = attributes.get_var_attributes();

    REQUIRE(variables.size() == 4);
    for (unsigned int index : {0, 1, 2, 3})
      {
        REQUIRE(variables.at(index).field_solve_type == fieldSolveType::EXPLICIT);
      }
  }

  /**
   * Two explicit equations that each have their own auxiliary variable (e.g., 2
   * Cahn-Hilliard equations).
   */
  SECTION("Two explicit two auxiliary")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public variableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      loadVariableAttributes() override
      {
        for (unsigned int index : {0, 1})
          {
            set_variable_name(index, "n" + std::to_string(index));
            set_variable_type(index, SCALAR);
            set_variable_equation_type(index, EXPLICIT_TIME_DEPENDENT);
            set_dependencies_value_term_RHS(index, "n" + std::to_string(index));
            set_dependencies_gradient_term_RHS(index,
                                               "grad(xi" + std::to_string(index) + ")");
          }
        for (unsigned int index : {2, 3})
          {
            set_variable_name(index, "xi" + std::to_string(index - 2));
            set_variable_type(index, SCALAR);
            set_variable_equation_type(index, AUXILIARY);
            set_dependencies_value_term_RHS(index, "n" + std::to_string(index - 2));
            set_dependencies_gradient_term_RHS(index,
                                               "grad(n" + std::to_string(index - 2) +
                                                 ")");
          }
      }
    };

    testVariableAttributeLoader attributes;
    attributes.init_variable_attributes();
    AttributesList variables = attributes.get_var_attributes();

    REQUIRE(variables.size() == 4);
    for (unsigned int index : {0, 1})
      {
        REQUIRE(variables.at(index).field_solve_type == fieldSolveType::EXPLICIT);
      }
    for (unsigned int index : {2, 3})
      {
        REQUIRE(variables.at(index).field_solve_type ==
                fieldSolveType::NONEXPLICIT_AUXILIARY);
      }
  }

  /**
   * Two explicit equations that dependent on two different time-independent quantities.
   * The second time-independent is coupled to the value of the other. In this case we
   * have two Allen-Cahn equations with two Poisson equations. Note that the second
   * Poisson equation has a forcing term equal to the first (e.g, Δu=f).
   */
  SECTION("Two explicit two time-independent")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public variableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      loadVariableAttributes() override
      {
        for (unsigned int index : {0, 1})
          {
            set_variable_name(index, "n" + std::to_string(index));
            set_variable_type(index, VECTOR);
            set_variable_equation_type(index, EXPLICIT_TIME_DEPENDENT);
            set_dependencies_value_term_RHS(index,
                                            "grad(phi), eta, n" + std::to_string(index));
            set_dependencies_gradient_term_RHS(index,
                                               "grad(n" + std::to_string(index) + ")");
          }

        set_variable_name(2, "phi");
        set_variable_type(2, SCALAR);
        set_variable_equation_type(2, TIME_INDEPENDENT);

        set_dependencies_value_term_LHS(2, "");
        set_dependencies_gradient_term_LHS(2, "grad(change(phi))");
        set_dependencies_value_term_RHS(2, "");
        set_dependencies_gradient_term_RHS(2, "grad(phi)");

        set_variable_name(3, "eta");
        set_variable_type(3, SCALAR);
        set_variable_equation_type(3, TIME_INDEPENDENT);

        set_dependencies_value_term_LHS(3, "");
        set_dependencies_gradient_term_LHS(3, "grad(change(eta))");
        set_dependencies_value_term_RHS(3, "grad(phi)");
        set_dependencies_gradient_term_RHS(3, "grad(eta)");
      }
    };

    testVariableAttributeLoader attributes;
    attributes.init_variable_attributes();
    AttributesList variables = attributes.get_var_attributes();

    REQUIRE(variables.size() == 4);
    for (unsigned int index : {0, 1})
      {
        REQUIRE(variables.at(index).field_solve_type == fieldSolveType::EXPLICIT);
      }
    for (unsigned int index : {2, 3})
      {
        REQUIRE(variables.at(index).field_solve_type ==
                fieldSolveType::NONEXPLICIT_LINEAR);
      }
  }
}