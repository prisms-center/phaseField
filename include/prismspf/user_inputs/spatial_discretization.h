// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/grid_refiner_criterion.h>
#include <prismspf/core/types.h>

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
  template <typename MeshType>
  void
  mark_periodic(const MeshType &mesh)
  {
    if constexpr (std::is_same_v<MeshType, Triangulation>)
      {
        for (const auto &val : periodicity_set)
          {
            dealii::GridTools::collect_periodic_faces(mesh,
                                                      std::get<0>(val),
                                                      std::get<1>(val),
                                                      std::get<2>(val),
                                                      triangulation_periodiciy_vector);
          }
      }
    else if constexpr (std::is_same_v<MeshType, dealii::DoFHandler<dim>>)
      {
        for (const auto &val : periodicity_set)
          {
            dealii::GridTools::collect_periodic_faces(mesh,
                                                      std::get<0>(val),
                                                      std::get<1>(val),
                                                      std::get<2>(val),
                                                      dof_handler_periodiciy_vector);
          }
      }
  };

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
    triangulation_periodiciy_vector;

  /**
   * @brief DoFHandler periodicity vector.
   */
  std::vector<
    dealii::GridTools::PeriodicFacePair<typename dealii::DoFHandler<dim>::cell_iterator>>
    dof_handler_periodiciy_vector;
};

/**
 * @brief Class for rectangular mesh parameters.
 */
template <unsigned int dim>
struct RectangularMesh
{
  /**
   * @brief Constructor.
   */
  RectangularMesh() = default;

  /**
   * @brief Constructor.
   */
  RectangularMesh(dealii::Tensor<1, dim, double> _size,
                  dealii::Tensor<1, dim, double> _lower_bound,
                  std::array<unsigned int, dim>  _subdivisions)
    : size(_size)
    , lower_bound(_lower_bound)
    , subdivisions(_subdivisions) {};

  /**
   * @brief Generate the mesh.
   */
  void
  generate_mesh(Triangulation<dim> &triangulation) const
  {
    validate();
    dealii::GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                      subdivisions,
                                                      dealii::Point<dim>(lower_bound),
                                                      dealii::Point<dim>(size));
    mark_boundaries(triangulation);
    mark_periodic(triangulation);
  };

  /**
   * @brief Validate
   */
  void
  validate() const
  {}

  /**
   * @brief Domain extents in each cartesian direction.
   */
  dealii::Tensor<1, dim, double> size;

  /**
   * @brief Domain extents in each cartesian direction.
   */
  dealii::Tensor<1, dim, double> lower_bound;

  /**
   * @brief Mesh subdivisions in each cartesian direction.
   */
  std::vector<unsigned int> subdivisions = std::vector<unsigned int>(dim, 1);

  /**
   * @brief Which directions have periodic conditions
   */
  std::set<unsigned int> periodic_directions;

  /**
   * @brief Mark the boundaries of the mesh.
   */
  void
  mark_boundaries(Triangulation<dim> &triangulation) const
  {
    // Loop through the cells
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        // Mark the faces (faces_per_cell = 2*dim)
        for (unsigned int face_number = 0;
             face_number < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face_number)
          {
            // Direction for quad and hex cells
            unsigned int direction = face_number / 2;

            // Mark the boundary id for x=0, y=0, z=0 and x=max, y=max, z=max
            if (std::fabs(cell->face(face_number)->center()(direction)) <
                  Defaults::mesh_tolerance ||
                std::fabs(cell->face(face_number)->center()(direction) -
                          size[direction]) < Defaults::mesh_tolerance)
              {
                cell->face(face_number)->set_boundary_id(face_number);
              }
          }
      }
  }

  /**
   * @brief Mark the periodic faces of the mesh.
   */
  void
  mark_periodic(Triangulation<dim> &triangulation) const
  {
    // Create a vector of matched pairs that we fill and enforce upon the
    // constraints
    std::vector<
      dealii::GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
      periodicity_vector;
    collect_periodic_faces(triangulation, periodicity_vector);

    // Add periodicity
    triangulation.add_periodicity(periodicity_vector);
  }

  /**
   * @brief Get the periodic face pairs
   */
  template <typename MeshType> // triangulation or dofhandler
  void
  collect_periodic_faces(
    const MeshType &triangulation,
    std::vector<dealii::GridTools::PeriodicFacePair<typename MeshType::cell_iterator>>
      &periodicity_vector) const
  {
    for (unsigned int direction : periodic_directions)
      {
        // Collect the matched pairs on the coarsest level of the mesh
        unsigned int boundary_id = direction * 2;
        dealii::GridTools::collect_periodic_faces(triangulation,
                                                  boundary_id,
                                                  boundary_id + 1,
                                                  direction,
                                                  periodicity_vector);
      }
  }
};

/**
 * @brief Class for spherical mesh parameters.
 */
template <unsigned int dim>
struct SphericalMesh
{
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
  generate_mesh(Triangulation<dim> &triangulation) const
  {
    validate();
    // TODO: can we just generate a 1d mesh using
    // dealii::GridGenerator::subdivided_hyper_rectangle instead of throwing an exception?
    AssertThrow(dim != 1, dealii::ExcMessage("Spherical mesh not valid in 1D"));
    dealii::GridGenerator::hyper_ball(triangulation, dealii::Point<dim>(), radius);
    mark_boundaries(triangulation);
  };

  /**
   * @brief Validate
   */
  void
  validate() const
  {}

  /**
   * @brief Radius of the spherical domain.
   */
  double radius = 1.0;

private:
  /**
   * @brief Mark the boundaries of the mesh.
   */
  void
  mark_boundaries([[maybe_unused]] Triangulation<dim> &triangulation) const
  {
    // TODO mark all as 0, or come up with something complicated
  }
};

/**
 * @brief Struct that holds spatial discretization parameters.
 */
template <unsigned int dim>
struct SpatialDiscretization
{
  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  /**
   * @brief Validate
   */
  void
  validate()
  { // Check that AMR is not enabled for 1D
    AssertThrow(
      (!has_adaptivity || dim != 1),
      dealii::ExcMessage(
        "Adaptive meshing for the matrix-free method is not currently supported."));

    // Some check if AMR is enabled
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
    // Validate whichever mesh generator is being used
    if (type == TriangulationType::Rectangular)
      {
        rectangular_mesh.validate();
      }
    else if (type == TriangulationType::Spherical)
      {
        spherical_mesh.validate();
      }
    else if (type == TriangulationType::Custom)
      {
        return;
      }
    AssertThrow(false, UnreachableCode("Invalid TriangulationType"));
  }

  /**
   * @brief Whether the provided increment is a valid grid refinement step.
   */
  [[nodiscard]] bool
  should_refine_mesh(unsigned int increment) const
  {
    return increment % remeshing_period == 0;
  }

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler,
                     unsigned int              max_criteria = 5) const;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler,
                    unsigned int              max_criteria = 5);

  // Triangulation type
  TriangulationType type = TriangulationType::Rectangular;

  // Rectangular mesh parameters
  RectangularMesh<dim> rectangular_mesh;

  // Spherical mesh parameters
  SphericalMesh<dim> spherical_mesh;

  // Global refinement of mesh
  unsigned int global_refinement = 0;

  // Element polynomial degree
  unsigned int degree = 1;

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

template <unsigned int dim>
inline void
SpatialDiscretization<dim>::print_parameter_summary() const
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Spatial Discretization\n"
    << "================================================\n";

  if (type == TriangulationType::Spherical)
    {
    }
  else if (type == TriangulationType::Rectangular)
    {
    }
  else if (type == TriangulationType::Custom)
    {
    }

  ConditionalOStreams::pout_summary()
    << "Global refinement: " << global_refinement << "\n"
    << "Degree: " << degree << "\n"
    << "Adaptivity enabled: " << bool_to_string(has_adaptivity) << "\n"
    << "Max refinement: " << max_refinement << "\n"
    << "Min refinement: " << min_refinement << "\n"
    << "Remeshing period: " << remeshing_period << "\n";

  if (!refinement_criteria.empty())
    {
      ConditionalOStreams::pout_summary() << "Refinement criteria:\n";
      for (const auto &[field_name, criterion] : refinement_criteria)
        {
          ConditionalOStreams::pout_summary()
            << "  Field name: " << field_name << "\n"
            << "  Criterion type: " << criterion.criterion_string() << "\n"
            << "  Value lower bound: " << criterion.value_lower_bound << "\n"
            << "  Value upper bound: " << criterion.value_upper_bound << "\n"
            << "  Gradient lower bound: " << criterion.gradient_lower_bound << "\n\n";
        }
    }

  ConditionalOStreams::pout_summary() << "\n" << std::flush;
}

template <unsigned int dim>
inline void
SpatialDiscretization<dim>::declare_parameters(
  dealii::ParameterHandler &parameter_handler,
  unsigned int              max_criteria) const
{
  parameter_handler.declare_entry("global refinement",
                                  "0",
                                  dealii::Patterns::Integer(0, INT_MAX),
                                  "The number of initial refinements of the coarse mesh.",
                                  true);

  parameter_handler.enter_subsection("Rectangular mesh");
  {
    parameter_handler.declare_entry("x size",
                                    "0.0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The size of the domain in the x direction.");
    parameter_handler.declare_entry("y size",
                                    "0.0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The size of the domain in the y direction.");
    parameter_handler.declare_entry("z size",
                                    "0.0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The size of the domain in the z direction.");
    parameter_handler.declare_entry("x lower bound",
                                    "0.0",
                                    dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                    "The lower bound of the domain in the x direction.");
    parameter_handler.declare_entry("y lower bound",
                                    "0.0",
                                    dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                    "The lower bound of the domain in the y direction.");
    parameter_handler.declare_entry("z lower bound",
                                    "0.0",
                                    dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                    "The lower bound of the domain in the z direction.");
    parameter_handler.declare_entry(
      "x subdivisions",
      "1",
      dealii::Patterns::Integer(1, INT_MAX),
      "The number of mesh subdivisions in the x direction.");
    parameter_handler.declare_entry(
      "y subdivisions",
      "1",
      dealii::Patterns::Integer(1, INT_MAX),
      "The number of mesh subdivisions in the y direction.");
    parameter_handler.declare_entry(
      "z subdivisions",
      "1",
      dealii::Patterns::Integer(1, INT_MAX),
      "The number of mesh subdivisions in the z direction.");

    parameter_handler.declare_entry("periodic x",
                                    "false",
                                    dealii::Patterns::Bool(),
                                    "Whether to have periodicity in the x-direction.");
    parameter_handler.declare_entry("periodic y",
                                    "false",
                                    dealii::Patterns::Bool(),
                                    "Whether to have periodicity in the y-direction.");
    parameter_handler.declare_entry("periodic z",
                                    "false",
                                    dealii::Patterns::Bool(),
                                    "Whether to have periodicity in the z-direction.");
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Spherical mesh");
  {
    parameter_handler.declare_entry("radius",
                                    "0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The radius of the domain.");
  }
  parameter_handler.leave_subsection();

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
          dealii::Patterns::Anything(),
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
}

template <unsigned int dim>
inline void
SpatialDiscretization<dim>::assign_parameters(dealii::ParameterHandler &parameter_handler,
                                              unsigned int              max_criteria)
{
  // Rectangular mesh
  parameter_handler.enter_subsection("Rectangular mesh");
  {
    static const std::vector<std::string> axis_labels = {"x", "y", "z"};
    RectangularMesh<dim>                 &rect        = rectangular_mesh;
    for (unsigned int i = 0; i < dim; ++i)
      {
        rect.size[i] = parameter_handler.get_double(axis_labels[i] + " size");
        rect.lower_bound[i] =
          parameter_handler.get_double(axis_labels[i] + " lower bound");
        rect.subdivisions[i] = static_cast<unsigned int>(
          parameter_handler.get_integer(axis_labels[i] + " subdivisions"));
        if (parameter_handler.get_bool("periodic " + axis_labels[i]))
          {
            rect.periodic_directions.insert(i);
          }
      }
  }
  parameter_handler.leave_subsection();
  // Spherical mesh
  parameter_handler.enter_subsection("Spherical mesh");
  {
    spherical_mesh.radius = parameter_handler.get_double("radius");
  }
  parameter_handler.leave_subsection();

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
        // todo: make case insensitive
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
}

PRISMS_PF_END_NAMESPACE
