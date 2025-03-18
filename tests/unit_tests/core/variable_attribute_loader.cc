// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/config.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>

#include "catch.hpp"

#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * This unit tests looks at variable_attributes.h and variable_attribute_loader.h and how
 * the parser interpretes invalid user dependencies.
 */
TEST_CASE("Invalid dependencies")
{
  SECTION("Dependencies on a constant explicit variable")
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
        set_variable_equation_type(0, CONSTANT);

        set_dependencies_value_term_RHS(0, "phi");
        set_dependencies_gradient_term_RHS(0, "");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("LHS dependencies on an explicit variable")
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
        set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

        set_dependencies_value_term_LHS(0, "phi");
        set_dependencies_gradient_term_LHS(0, "grad(phi)");
        set_dependencies_value_term_RHS(0, "phi");
        set_dependencies_gradient_term_RHS(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Empty variable name")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public variableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      loadVariableAttributes() override
      {
        set_variable_name(0, "");
        set_variable_type(0, SCALAR);
        set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Invalid field type")
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
        set_variable_type(0, UNDEFINED_FIELD);
        set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

        set_dependencies_value_term_RHS(0, "phi");
        set_dependencies_gradient_term_RHS(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Invalid pde type")
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
        set_variable_equation_type(0, UNDEFINED_PDE);

        set_dependencies_value_term_RHS(0, "phi");
        set_dependencies_gradient_term_RHS(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Invalid dependency name")
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
        set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

        set_dependencies_value_term_RHS(0, "fi");
        set_dependencies_gradient_term_RHS(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("TIME_INDEPENDENT postprocess variable")
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
        set_is_postprocessed_field(0, true);

        set_dependencies_value_term_LHS(0, "phi");
        set_dependencies_gradient_term_LHS(0, "grad(phi)");
        set_dependencies_value_term_RHS(0, "phi");
        set_dependencies_gradient_term_RHS(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("AUXILIARY postprocess variable")
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
        set_variable_equation_type(0, AUXILIARY);
        set_is_postprocessed_field(0, true);

        set_dependencies_value_term_RHS(0, "phi");
        set_dependencies_gradient_term_RHS(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("CONSTANT postprocess variable")
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
        set_variable_equation_type(0, CONSTANT);
        set_is_postprocessed_field(0, true);
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("IMPLICIT_TIME_DEPENDENT postprocess variable")
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
        set_variable_equation_type(0, IMPLICIT_TIME_DEPENDENT);
        set_is_postprocessed_field(0, true);

        set_dependencies_value_term_LHS(0, "phi");
        set_dependencies_gradient_term_LHS(0, "grad(phi)");
        set_dependencies_value_term_RHS(0, "phi");
        set_dependencies_gradient_term_RHS(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Sequential old solution dependencies")
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
        set_variable_equation_type(0, IMPLICIT_TIME_DEPENDENT);

        set_dependencies_value_term_LHS(0, "phi");
        set_dependencies_gradient_term_LHS(0, "grad(phi)");
        set_dependencies_value_term_RHS(0, "phi, old_4(phi)");
        set_dependencies_gradient_term_RHS(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Postprocess dependencies for main variable")
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
        set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

        set_dependencies_value_term_RHS(0, "phi, free_energy");
        set_dependencies_gradient_term_RHS(0, "grad(phi)");

        set_variable_name(1, "free_energy");
        set_variable_type(1, SCALAR);
        set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);
        set_is_postprocessed_field(1, true);

        set_dependencies_value_term_RHS(1, "phi");
        set_dependencies_gradient_term_RHS(1, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Postprocess dependencies for postprocess variable")
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
        set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

        set_dependencies_value_term_RHS(0, "phi");
        set_dependencies_gradient_term_RHS(0, "grad(phi)");

        set_variable_name(1, "free_energy");
        set_variable_type(1, SCALAR);
        set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);
        set_is_postprocessed_field(1, true);

        set_dependencies_value_term_RHS(1, "phi");
        set_dependencies_gradient_term_RHS(1, "grad(phi)");

        set_variable_name(2, "unavailable_energy");
        set_variable_type(2, SCALAR);
        set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);
        set_is_postprocessed_field(2, true);

        set_dependencies_value_term_RHS(2, "free_energy");
        set_dependencies_gradient_term_RHS(2, "grad(free_energy)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }
}

PRISMS_PF_END_NAMESPACE