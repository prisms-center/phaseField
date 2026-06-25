// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/affine_constraints.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/grid_refiner_criterion.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Internal enum for various triangulation types.
 */
enum TriangulationType : std::uint8_t
{
  Rectangular,
  Spherical,
  Custom
};

/**
 * @brief Triangulation base class.
 */
template <unsigned int dim>
struct Mesh
{
  using Triangulation =
    std::conditional_t<dim == 1,
                       dealii::Triangulation<dim>,
                       dealii::parallel::distributed::Triangulation<dim>>;

  /**
   * @brief Constructor.
   */
  Mesh() = default;

  /**
   * @brief Generate the mesh.
   */
  virtual void
  generate_mesh(Triangulation &triangulation) const = 0;

  /**
   * @brief Mark the boundaries of the mesh.
   */
  virtual void
  mark_boundaries(Triangulation &triangulation) const;

  /**
   * @brief Mark the periodic faces of the mesh.
   */
  void
  mark_periodic(Triangulation &tria)
  {
    // Loop over provided periodicity set and add to the periodicity vector
    for (const auto &[id_1, id_2, direction] : periodicity_set)
      {
        dealii::GridTools::collect_periodic_faces(tria,
                                                  id_1,
                                                  id_2,
                                                  direction,
                                                  triangulation_periodicity_vector);
      }
    // Pass periodicity vector to triangulation
    tria.add_periodicity(triangulation_periodicity_vector);
  };

  /**
   * @brief Mark the periodic faces of the mesh.
   */
  template <typename number>
  void
  mark_periodic(dealii::DoFHandler<dim>           &dof_handler,
                dealii::AffineConstraints<number> &constraints)
  {
    // Loop over provided periodicity set and add to the periodicity vector
    for (const auto &[id_1, id_2, direction] : periodicity_set)
      {
        dealii::GridTools::collect_periodic_faces(dof_handler,
                                                  id_1,
                                                  id_2,
                                                  direction,
                                                  dof_handler);
      }

    // Pass periodicity vector to constraints
    // NOTE: We may want to add a mask here. If we're dealing with vector fields, this
    // will enforce periodicity constraints for all components on a given face.
    dealii::DoFTools::make_periodicity_constraints(dof_handler_periodicity_vector,
                                                   constraints);
  };

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  virtual void
  declare_parameters(dealii::ParameterHandler &parameter_handler) const = 0;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  virtual void
  assign_parameters(dealii::ParameterHandler &parameter_handler) = 0;

  /**
   * @brief Validate.
   */
  virtual void
  validate() const = 0;

  /**
   * @brief Periodicity set.
   *
   * This set contains all the information about periodicity for the system. Each entry
   * contain tuple with three entries: the first boundary id, the second boundary id, and
   * the direction.
   *
   * For example, (0,1,0) would link the 0th and 1st boundaries with periodicity in the
   * x-direction.
   */
  std::set<std::tuple<unsigned int, unsigned int, unsigned int>> periodicity_set;

  /**
   * @brief Triangulation periodicity vector.
   */
  std::vector<dealii::GridTools::PeriodicFacePair<typename Triangulation::cell_iterator>>
    triangulation_periodicity_vector;

  /**
   * @brief DoFHandler periodicity vector.
   */
  std::vector<
    dealii::GridTools::PeriodicFacePair<typename dealii::DoFHandler<dim>::cell_iterator>>
    dof_handler_periodicity_vector;
};

/**
 * @brief Class for rectangular mesh parameters.
 */
template <unsigned int dim>
struct RectangularMesh : public Mesh<dim>
{
  using Triangulation = typename Mesh<dim>::Triangulation;

  /**
   * @brief Constructor.
   */
  RectangularMesh() = default;

  /**
   * @brief Constructor.
   */
  RectangularMesh(dealii::Tensor<1, dim, double> _upper_bound,
                  dealii::Tensor<1, dim, double> _lower_bound,
                  std::array<unsigned int, dim>  _subdivisions)
    : upper_bound(_upper_bound)
    , lower_bound(_lower_bound)
    , subdivisions(_subdivisions) {};

  /**
   * @brief Generate the mesh.
   */
  void
  generate_mesh(Triangulation &triangulation) const override
  {
    validate();
    dealii::GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                      subdivisions,
                                                      dealii::Point<dim>(lower_bound),
                                                      dealii::Point<dim>(upper_bound),
                                                      true);
  };

  /**
   * @brief Mark the boundaries of the mesh.
   */
  void
  mark_boundaries(Triangulation &triangulation) const override {
    // The colorize option above does this for us.
    //
    // Here are the mappings for reference,
    //
    // x=0 -> 0
    // x=max -> 1
    // y=0 -> 2
    // y=max -> 3
    // z=0 -> 4
    // z=max -> 5
  };

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler) const override
  {
    parameter_handler.enter_subsection("Rectangular mesh");
    {
      for (const auto &dir : {"x", "y", "z"})
        {
          parameter_handler.declare_entry(std::string(dir) + " upper bound",
                                          "0.0",
                                          dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                          "The upper bound of the domain in the " +
                                            std::string(dir) + "-direction.");
          parameter_handler.declare_alias(std::string(dir) + " upper bound",
                                          std::string(dir) + " size",
                                          true);
          parameter_handler.declare_alias(std::string(dir) + " upper bound",
                                          "upper bound" + std::string(dir));

          parameter_handler.declare_entry(std::string(dir) + " lower bound",
                                          "0.0",
                                          dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                          "The lower bound of the domain in the " +
                                            std::string(dir) + "-direction.");
          parameter_handler.declare_alias(std::string(dir) + " lower bound",
                                          "lower bound" + std::string(dir));

          parameter_handler.declare_entry(std::string(dir) + " subdivisions",
                                          "1",
                                          dealii::Patterns::Integer(1, INT_MAX),
                                          "The number of mesh subdivisions in the " +
                                            std::string(dir) + "-direction.");
          parameter_handler.declare_alias(std::string(dir) + " subdivisions",
                                          "subdivisions" + std::string(dir));

          parameter_handler.declare_entry(std::string(dir) + " periodic",
                                          "false",
                                          dealii::Patterns::Bool(),
                                          "Whether to have periodicity in the " +
                                            std::string(dir) + "-direction.");
          parameter_handler.declare_alias(std::string(dir) + " periodic",
                                          "periodic" + std::string(dir));
        }
    }
    parameter_handler.leave_subsection();
  };

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler) override
  {
    parameter_handler.enter_subsection("Rectangular mesh");
    {
      static const std::array<std::string_view, 3> axis_labels {"x", "y", "z"};

      for (unsigned int i = 0; i < dim; ++i)
        {
          const std::string axis {axis_labels.at(i)};

          upper_bound[i] = parameter_handler.get_double(axis + " upper bound");

          lower_bound[i] = parameter_handler.get_double(axis + " lower bound");

          subdivisions[i] = static_cast<unsigned int>(
            parameter_handler.get_integer(axis + " subdivisions"));

          if (parameter_handler.get_bool(axis + " periodic"))
            {
              this->periodicity_set.insert({2 * i, 2 * i + 1, i});
            }
        }
    }
    parameter_handler.leave_subsection();
  };

  /**
   * @brief Validate
   */
  void
  validate() const override
  {
    AssertThrow(
      (upper_bound - lower_bound).norm() != 0.0,
      dealii::ExcMessage(
        "Upper and lower bound for the mesh are the same point (total size is 0)."));
    // These next two asserts should be caught earlier if users are using the parameters
    // file. This is mostly for users that are using the lower level structures.
    Assert(
      subdivisions.size() == dim,
      dealii::ExcMessage(
        "Subdivisions vector size is not equal to the number of cartesian directions."));
    for (const auto &val : subdivisions)
      {
        Assert(val > 0, dealii::ExcMessage("Subdivisions must be greater than 0."));
      }
  }

  /**
   * @brief Upper bound point.
   */
  dealii::Tensor<1, dim, double> upper_bound;

  /**
   * @brief Lower bound point.
   */
  dealii::Tensor<1, dim, double> lower_bound;

  /**
   * @brief Mesh subdivisions in each cartesian direction.
   */
  std::vector<unsigned int> subdivisions = std::vector<unsigned int>(dim, 1);
};

/**
 * @brief Class for spherical mesh parameters.
 */
template <unsigned int dim>
struct SphericalMesh : public Mesh<dim>
{
  using Triangulation = typename Mesh<dim>::Triangulation;

  /**
   * @brief Constructor.
   */
  SphericalMesh() = default;

  /**
   * @brief Constructor.
   */
  explicit SphericalMesh(double _radius)
    : radius(_radius)
  {}

  /**
   * @brief Generate the mesh.
   */
  void
  generate_mesh(Triangulation &triangulation) const
  {
    validate();
    // It doesn't make sense to use spherical meshes in 1D. Users should just switch to
    // rectangular.
    //
    // We could just rectangular here, but then we would run into issues where the coarse
    // mesh might be different between the hyper ball and the subdivided rectangle. Better
    // off to have users change it in their parameter file. It's only a few lines anyway.
    AssertThrow(dim != 1, dealii::ExcMessage("Spherical mesh not valid in 1D"));
    dealii::GridGenerator::hyper_ball(triangulation, dealii::Point<dim>(), radius);
  };

  /**
   * @brief Mark the boundaries of the mesh.
   */
  void
  mark_boundaries(Triangulation &triangulation) const override
  {
    // There's only 1 boundary on a sphere
  }

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler) const override
  {
    parameter_handler.enter_subsection("Spherical mesh");
    {
      parameter_handler.declare_entry("radius",
                                      "0",
                                      dealii::Patterns::Double(0.0, DBL_MAX),
                                      "The radius of the domain.");
    }
    parameter_handler.leave_subsection();
  };

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler) override
  {
    // Spherical mesh
    parameter_handler.enter_subsection("Spherical mesh");
    {
      radius = parameter_handler.get_double("radius");
    }
    parameter_handler.leave_subsection();
  };

  /**
   * @brief Validate
   */
  void
  validate() const override
  {
    AssertThrow(radius > 0.0, dealii::ExcMessage("Radius must be greater than zero."));
  }

  /**
   * @brief Radius of the spherical domain.
   */
  double radius = 0.0;
};

/**
 * @brief Struct that holds spatial discretization parameters.
 */
template <unsigned int dim>
struct SpatialDiscretization : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  virtual void
  predeclare(dealii::ParameterHandler &parameter_handler) const
  {
    parameter_handler.declare_entry("mesh type",
                                    "rectangular",
                                    dealii::Patterns::Selection(
                                      "rectangular|spherical|custom"),
                                    "The type of mesh to use.",
                                    true);
  };

  /**
   * @brief Assign the parameters from file.
   */
  virtual void
  preassign(dealii::ParameterHandler &parameter_handler)
  {
    const std::string mesh_type = parameter_handler.get("mesh type");

    if (mesh_type == "rectangular")
      {
        mesh = std::make_unique<RectangularMesh<dim>>();
      }
    else if (mesh_type == "spherical")
      {
        mesh = std::make_unique<SphericalMesh<dim>>();
      }
  };

  /**
   * @brief Declare the parameters to be read from file.
   */
  virtual void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int              max_criteria = Numbers::max_subsections) const
  {
    mesh->declare_parameters(parameter_handler);

    parameter_handler.declare_entry(
      "global refinement",
      "0",
      dealii::Patterns::Integer(0, INT_MAX),
      "The number of initial refinements of the coarse mesh.",
      true);

    parameter_handler.declare_entry("mesh adaptivity",
                                    "false",
                                    dealii::Patterns::Bool(),
                                    "Whether to enable mesh adaptivity.");
    parameter_handler.declare_entry("max refinement",
                                    "0",
                                    dealii::Patterns::Integer(0, INT_MAX),
                                    "The maximum level of refinement.");
    parameter_handler.declare_entry("min refinement",
                                    "0",
                                    dealii::Patterns::Integer(0, INT_MAX),
                                    "The minimum level of refinement.");
    parameter_handler.declare_entry(
      "remeshing period",
      "2147483647",
      dealii::Patterns::Integer(1, INT_MAX),
      "The number of time steps between mesh refinement operations.");

    for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
      {
        std::string subsection_text =
          "refinement criterion: " + std::to_string(criterion_id);
        parameter_handler.enter_subsection(subsection_text);
        {
          parameter_handler.declare_entry(
            "variables",
            "",
            dealii::Patterns::List(dealii::Patterns::Anything(), 0, INT_MAX, ","),
            "The names of the fields that will use this refinement criterion.");
          parameter_handler.declare_entry(
            "type",
            "none",
            dealii::Patterns::Selection("none|value|gradient|value_and_gradient"),
            "The type of criterion used to determine if a cell should be "
            "refined. The options are none, value, gradient, value_and_gradient.");
          parameter_handler.declare_entry(
            "value lower bound",
            "0.0",
            dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
            "The lower bound for the window determining where the mesh should be "
            "refined.");
          parameter_handler.declare_entry(
            "value upper bound",
            "0.0",
            dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
            "The upper bound for the window determining where the mesh should be "
            "refined.");
          parameter_handler.declare_entry("gradient magnitude lower bound",
                                          "2147483647",
                                          dealii::Patterns::Double(0.0, DBL_MAX),
                                          "The magnitude of the gradient above "
                                          "which the mesh should be refined.");
          parameter_handler.declare_alias("gradient magnitude lower bound",
                                          "gradient lower bound");
        }
        parameter_handler.leave_subsection();
      }
  };

  /**
   * @brief Assign the parameters from file.
   */
  virtual void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections)
  {
    mesh->assign_parameters(parameter_handler);

    global_refinement =
      (static_cast<unsigned int>(parameter_handler.get_integer("global refinement")));

    has_adaptivity = (parameter_handler.get_bool("mesh adaptivity"));

    remeshing_period =
      (static_cast<unsigned int>(parameter_handler.get_integer("remeshing period")));

    max_refinement =
      (static_cast<unsigned int>(parameter_handler.get_integer("max refinement")));
    min_refinement =
      (static_cast<unsigned int>(parameter_handler.get_integer("min refinement")));

    for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
      {
        std::string subsection_text =
          "refinement criterion: " + std::to_string(criterion_id);
        parameter_handler.enter_subsection(subsection_text);
        {
          std::vector<std::string> field_names =
            dealii::Utilities::split_string_list(parameter_handler.get("variables"));
          static const std::map<std::string, RefinementFlags> crit_map {
            {"",                   RefinementFlags::Nothing                          },
            {"none",               RefinementFlags::Nothing                          },
            {"value",              RefinementFlags::Value                            },
            {"gradient",           RefinementFlags::Gradient                         },
            {"value_and_gradient", RefinementFlags::Value | RefinementFlags::Gradient}
          };
          // TODO: make case insensitive
          RefinementCriterion criterion(crit_map.at(parameter_handler.get("type")),
                                        parameter_handler.get_double("value lower bound"),
                                        parameter_handler.get_double("value upper bound"),
                                        parameter_handler.get_double(
                                          "gradient magnitude lower bound"));
          for (const auto &field_name : field_names)
            {
              refinement_criteria[field_name] = criterion;
            }
        }
        parameter_handler.leave_subsection();
      }
  };

  /**
   * @brief Validate.
   */
  void
  validate([[maybe_unused]] const std::vector<FieldAttributes> &field_attributes,
           [[maybe_unused]] const std::vector<SolveBlock> &solve_blocks) const override
  {
    mesh->validate();

    // Check that AMR is not enabled for 1D
    AssertThrow(
      (!has_adaptivity || dim != 1),
      dealii::ExcMessage(
        "Adaptive meshing for the matrix-free method is not currently supported."));

    // Some checks if AMR is enabled
    if (has_adaptivity)
      {
        // Check that the minimum and maximum refinement levels are valid
        AssertThrow((max_refinement >= global_refinement),
                    dealii::ExcMessage(
                      "The maximum refinement level must be greater than or equal to the "
                      "global refinement level."));
        AssertThrow((min_refinement <= global_refinement),
                    dealii::ExcMessage(
                      "The minimum refinement level must be less than or equal to the "
                      "global refinement level."));
        AssertThrow((max_refinement >= min_refinement),
                    dealii::ExcMessage(
                      "The maximum refinement level must be greater than or equal to the "
                      "minimum refinement level."));
      }
  }

  /**
   * @brief Whether the provided increment is a valid grid refinement step.
   */
  [[nodiscard]] bool
  should_refine_mesh(unsigned int increment) const
  {
    return increment % remeshing_period == 0;
  }

  // Mesh object
  std::unique_ptr<Mesh<dim>> mesh;

  // Global refinement of mesh
  unsigned int global_refinement = 0;

  // Whether adaptive meshing (AMR) is enabled
  bool has_adaptivity = false;

  // Maximum global refinement for AMR
  unsigned int max_refinement = 0;

  // Minimum global refinement for AMR
  unsigned int min_refinement = 0;

  // The number of steps between remeshing
  unsigned int remeshing_period = UINT_MAX;

  // The criteria used for remeshing
  std::map<std::string, RefinementCriterion> refinement_criteria;
};

PRISMS_PF_END_NAMESPACE
