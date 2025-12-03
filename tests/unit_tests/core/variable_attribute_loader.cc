// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

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
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(0, Constant);

        set_dependencies_value_term_rhs(0, "phi");
        set_dependencies_gradient_term_rhs(0, "");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("LHS dependencies on an explicit variable")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(0, ExplicitTimeDependent);

        set_dependencies_value_term_lhs(0, "phi");
        set_dependencies_gradient_term_lhs(0, "grad(phi)");
        set_dependencies_value_term_rhs(0, "phi");
        set_dependencies_gradient_term_rhs(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Empty variable name")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "");
        set_variable_type(0, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(0, ExplicitTimeDependent);
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Invalid field type")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Undefined);
        set_variable_equation_type(0, ExplicitTimeDependent);

        set_dependencies_value_term_rhs(0, "phi");
        set_dependencies_gradient_term_rhs(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Invalid pde type")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Vector);
        set_variable_equation_type(0, Numbers::invalid_pde_type);

        set_dependencies_value_term_rhs(0, "phi");
        set_dependencies_gradient_term_rhs(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Invalid dependency name")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(0, ExplicitTimeDependent);

        set_dependencies_value_term_rhs(0, "fi");
        set_dependencies_gradient_term_rhs(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("TimeIndependent postprocess variable")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(0, TimeIndependent);
        set_is_postprocessed_field(0, true);

        set_dependencies_value_term_lhs(0, "phi");
        set_dependencies_gradient_term_lhs(0, "grad(phi)");
        set_dependencies_value_term_rhs(0, "phi");
        set_dependencies_gradient_term_rhs(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Auxiliary postprocess variable")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(0, Auxiliary);
        set_is_postprocessed_field(0, true);

        set_dependencies_value_term_rhs(0, "phi");
        set_dependencies_gradient_term_rhs(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Constant postprocess variable")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(0, Constant);
        set_is_postprocessed_field(0, true);
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("ImplicitTimeDependent postprocess variable")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(0, ImplicitTimeDependent);
        set_is_postprocessed_field(0, true);

        set_dependencies_value_term_lhs(0, "phi");
        set_dependencies_gradient_term_lhs(0, "grad(phi)");
        set_dependencies_value_term_rhs(0, "phi");
        set_dependencies_gradient_term_rhs(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Sequential old solution dependencies")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(0, ImplicitTimeDependent);

        set_dependencies_value_term_lhs(0, "phi");
        set_dependencies_gradient_term_lhs(0, "grad(phi)");
        set_dependencies_value_term_rhs(0, "phi, old_4(phi)");
        set_dependencies_gradient_term_rhs(0, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Postprocess dependencies for main variable")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(0, ExplicitTimeDependent);

        set_dependencies_value_term_rhs(0, "phi, free_energy");
        set_dependencies_gradient_term_rhs(0, "grad(phi)");

        set_variable_name(1, "free_energy");
        set_variable_type(1, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(1, ExplicitTimeDependent);
        set_is_postprocessed_field(1, true);

        set_dependencies_value_term_rhs(1, "phi");
        set_dependencies_gradient_term_rhs(1, "grad(phi)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }

  SECTION("Postprocess dependencies for postprocess variable")
  {
    // Create test class for variable attribute loader
    class testVariableAttributeLoader : public VariableAttributeLoader
    {
    public:
      ~testVariableAttributeLoader() override = default;

      void
      load_variable_attributes() override
      {
        set_variable_name(0, "phi");
        set_variable_type(0, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(0, ExplicitTimeDependent);

        set_dependencies_value_term_rhs(0, "phi");
        set_dependencies_gradient_term_rhs(0, "grad(phi)");

        set_variable_name(1, "free_energy");
        set_variable_type(1, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(1, ExplicitTimeDependent);
        set_is_postprocessed_field(1, true);

        set_dependencies_value_term_rhs(1, "phi");
        set_dependencies_gradient_term_rhs(1, "grad(phi)");

        set_variable_name(2, "unavailable_energy");
        set_variable_type(2, FieldInfo::TensorRank::Scalar);
        set_variable_equation_type(2, ExplicitTimeDependent);
        set_is_postprocessed_field(2, true);

        set_dependencies_value_term_rhs(2, "free_energy");
        set_dependencies_gradient_term_rhs(2, "grad(free_energy)");
      }
    };

    testVariableAttributeLoader attributes;
    REQUIRE_THROWS(attributes.init_variable_attributes());
  }
}

PRISMS_PF_END_NAMESPACE
