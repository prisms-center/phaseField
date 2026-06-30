#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/solve_block.h>

#include <prismspf/user_inputs/parameter_base.h>
#include <prismspf/user_inputs/spatial_discretization.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
void
Mesh<dim>::mark_periodic(typename Mesh<dim>::Triangulation &triangulation)
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
}

template <unsigned int dim>
template <typename number>
void
Mesh<dim>::mark_periodic(dealii::DoFHandler<dim>           &dof_handler,
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
}

template <unsigned int dim>
RectangularMesh<dim>::RectangularMesh(dealii::Tensor<1, dim, double> _upper_bound,
                                      dealii::Tensor<1, dim, double> _lower_bound,
                                      std::array<unsigned int, dim>  _subdivisions)
  : upper_bound(_upper_bound)
  , lower_bound(_lower_bound)
  , subdivisions(_subdivisions)
{}

template <unsigned int dim>
void
RectangularMesh<dim>::generate_mesh(
  typename RectangularMesh<dim>::Triangulation &triangulation) const
{
  validate();
  dealii::GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                    subdivisions,
                                                    dealii::Point<dim>(lower_bound),
                                                    dealii::Point<dim>(upper_bound),
                                                    true);
}

template <unsigned int dim>
void
RectangularMesh<dim>::mark_boundaries(
  typename RectangularMesh<dim>::Triangulation &triangulation) const
{
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
}

template <unsigned int dim>
void
RectangularMesh<dim>::declare_parameters(
  dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.enter_subsection("Rectangular mesh");
  {
    for (const auto &dir : axis_labels)
      {
        const std::string axis {dir};

        parameter_handler.declare_entry(axis + " upper bound",
                                        "0.0",
                                        dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                        "The upper bound of the domain in the " + axis +
                                          "-direction.");
        parameter_handler.declare_alias(axis + " upper bound", axis + " size", true);
        parameter_handler.declare_alias(axis + " upper bound", "upper bound" + axis);

        parameter_handler.declare_entry(axis + " lower bound",
                                        "0.0",
                                        dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                        "The lower bound of the domain in the " + axis +
                                          "-direction.");
        parameter_handler.declare_alias(axis + " lower bound", "lower bound" + axis);

        parameter_handler.declare_entry(axis + " subdivisions",
                                        "1",
                                        dealii::Patterns::Integer(1, INT_MAX),
                                        "The number of mesh subdivisions in the " + axis +
                                          "-direction.");
        parameter_handler.declare_alias(axis + " subdivisions", "subdivisions" + axis);

        parameter_handler.declare_entry(axis + " periodic",
                                        "false",
                                        dealii::Patterns::Bool(),
                                        "Whether to have periodicity in the " + axis +
                                          "-direction.");
        parameter_handler.declare_alias(axis + " periodic", "periodic" + axis);
      }
  }
  parameter_handler.leave_subsection();
}

template <unsigned int dim>
void
RectangularMesh<dim>::assign_parameters(dealii::ParameterHandler &parameter_handler)
{
  parameter_handler.enter_subsection("Rectangular mesh");
  {
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
}

template <unsigned int dim>
void
RectangularMesh<dim>::validate() const
{
  AssertThrow(
    (upper_bound - lower_bound).norm() != 0.0,
    dealii::ExcMessage(
      "Upper and lower bound for the mesh are the same point (total size is 0)."));
  // These next two asserts should be caught earlier if users are using the parameters
  // file. This is mostly for users that are using the lower level structures.
  Assert(subdivisions.size() == dim,
         dealii::ExcMessage("Subdivisions vector size is not equal to the number of "
                            "cartesian directions."));
  for (const auto &val : subdivisions)
    {
      Assert(val > 0, dealii::ExcMessage("Subdivisions must be greater than 0."));
    }
}

template <unsigned int dim>
SphericalMesh<dim>::SphericalMesh(double _radius)
  : radius(_radius)
{}

template <unsigned int dim>
void
SphericalMesh<dim>::generate_mesh(
  typename SphericalMesh<dim>::Triangulation &triangulation) const
{
  validate();
  // It doesn't make sense to use spherical meshes in 1D. Users should just switch to
  // rectangular.
  //
  // We could just rectangular here, but then we would run into issues where the
  // coarse mesh might be different between the hyper ball and the subdivided
  // rectangle. Better off to have users change it in their parameter file. It's only
  // a few lines anyway.
  AssertThrow(dim != 1, dealii::ExcMessage("Spherical mesh not valid in 1D"));
  dealii::GridGenerator::hyper_ball(triangulation, dealii::Point<dim>(), radius);
}

template <unsigned int dim>
void
SphericalMesh<dim>::mark_boundaries(
  typename SphericalMesh<dim>::Triangulation &triangulation) const
{
  // There's only 1 boundary on a sphere
}

template <unsigned int dim>
void
SphericalMesh<dim>::declare_parameters(dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.enter_subsection("Spherical mesh");
  {
    parameter_handler.declare_entry("radius",
                                    "0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The radius of the domain.");
  }
  parameter_handler.leave_subsection();
}

template <unsigned int dim>
void
SphericalMesh<dim>::assign_parameters(dealii::ParameterHandler &parameter_handler)
{
  parameter_handler.enter_subsection("Spherical mesh");
  {
    radius = parameter_handler.get_double("radius");
  }
  parameter_handler.leave_subsection();
}

template <unsigned int dim>
void
SphericalMesh<dim>::validate() const
{
  AssertThrow(radius > 0.0, dealii::ExcMessage("Radius must be greater than zero."));
}

template <unsigned int dim>
void
SpatialDiscretization<dim>::predeclare(dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.declare_entry("mesh type",
                                  "rectangular",
                                  dealii::Patterns::Selection(
                                    "rectangular|spherical|custom"),
                                  "The type of mesh to use.",
                                  true);
}

template <unsigned int dim>
void
SpatialDiscretization<dim>::preassign(dealii::ParameterHandler &parameter_handler)
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
  else
    {
      // If custom, it's up to the user to define the class and create the pointer.
    }
}

template <unsigned int dim>
void
SpatialDiscretization<dim>::declare(dealii::ParameterHandler &parameter_handler,
                                    unsigned int              max_criteria) const
{
  mesh->declare_parameters(parameter_handler);

  parameter_handler.declare_entry("global refinement",
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
}

template <unsigned int dim>
void
SpatialDiscretization<dim>::assign(dealii::ParameterHandler &parameter_handler,
                                   unsigned int              max_criteria)
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
}

template <unsigned int dim>
void
SpatialDiscretization<dim>::validate(
  [[maybe_unused]] const std::vector<FieldAttributes> &field_attributes,
  [[maybe_unused]] const std::vector<SolveBlock>      &solve_blocks) const
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

template <unsigned int dim>
bool
SpatialDiscretization<dim>::should_refine_mesh(unsigned int increment) const
{
  return increment % remeshing_period == 0;
}

PRISMS_PF_END_NAMESPACE
