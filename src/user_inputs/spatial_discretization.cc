#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/solve_block.h>

#include <prismspf/user_inputs/parameter_base.h>
#include <prismspf/user_inputs/spatial_discretization.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
double
Mesh<dim>::distance(const dealii::Point<dim> &point_1,
                    const dealii::Point<dim> &point_2) const
{
  return point_1.distance(point_2);
}

template <unsigned int dim>
std::list<PeriodicPair<dim>>
Mesh<dim>::periodicity_set() const
{
  return {};
}

template <unsigned int dim>
void
Mesh<dim>::mark_periodic(typename Mesh<dim>::Triangulation &triangulation) const
{
  // Triangulation periodicity vector.
  std::vector<dealii::GridTools::PeriodicFacePair<typename Triangulation::cell_iterator>>
    triangulation_periodicity_vector;
  // Loop over provided periodicity set and add to the periodicity vector
  for (const auto &periodic_pair : this->periodicity_set())
    {
      dealii::GridTools::collect_periodic_faces(triangulation,
                                                periodic_pair.boundary_id_1,
                                                periodic_pair.boundary_id_2,
                                                periodic_pair.direction,
                                                triangulation_periodicity_vector,
                                                periodic_pair.translation_vector,
                                                periodic_pair.rotation_matrix);
    }
  // Pass periodicity vector to triangulation
  triangulation.add_periodicity(triangulation_periodicity_vector);
}

template <unsigned int dim>
template <typename number>
void
Mesh<dim>::mark_periodic(const dealii::DoFHandler<dim>     &dof_handler,
                         dealii::AffineConstraints<number> &constraints) const
{
  // DoFHandler periodicity vector.
  std::vector<
    dealii::GridTools::PeriodicFacePair<typename dealii::DoFHandler<dim>::cell_iterator>>
    dof_handler_periodicity_vector;
  // Loop over provided periodicity set and add to the periodicity vector
  for (const auto &periodic_pair : this->periodicity_set())
    {
      dealii::GridTools::collect_periodic_faces(dof_handler,
                                                periodic_pair.boundary_id_1,
                                                periodic_pair.boundary_id_2,
                                                periodic_pair.direction,
                                                dof_handler_periodicity_vector,
                                                periodic_pair.translation_vector,
                                                periodic_pair.rotation_matrix);
    }

  // Pass periodicity vector to constraints
  // NOTE: We may want to add a mask here. If we're dealing with vector fields, this
  // will enforce periodicity constraints for all components on a given face.
  dealii::DoFTools::make_periodicity_constraints<dim, dim, number>(
    dof_handler_periodicity_vector,
    constraints);
}

template <unsigned int dim>
RectangularMesh<dim>::RectangularMesh(dealii::Tensor<1, dim, double> _size,
                                      dealii::Tensor<1, dim, double> _lower_bound,
                                      std::vector<unsigned int>      _subdivisions)
  : size(_size)
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
                                                    dealii::Point<dim>(lower_bound +
                                                                       size),
                                                    true);
}

template <unsigned int dim>
void
RectangularMesh<dim>::mark_boundaries(
  typename RectangularMesh<dim>::Triangulation &triangulation) const
{
  // The colorize option in `dealii::GridGenerator::subdivided_hyper_rectangle` does this
  // for us.
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
std::list<PeriodicPair<dim>>
RectangularMesh<dim>::periodicity_set() const
{
  std::list<PeriodicPair<dim>> p_set;
  for (const auto &dir : periodic_directions)
    {
      PeriodicPair<dim> p_pair;
      p_pair.boundary_id_1 = 2 * dir;
      p_pair.boundary_id_2 = 2 * dir + 1;
      p_pair.direction     = dir;
      p_set.push_back(p_pair);
    }
  return p_set;
}

template <unsigned int dim>
double
RectangularMesh<dim>::distance(const dealii::Point<dim> &point_1,
                               const dealii::Point<dim> &point_2) const
{
  using std::sqrt;
  const std::list<PeriodicPair<dim>> &pair_set = periodicity_set();
  if (pair_set.empty())
    {
      return point_1.distance(point_2);
    }

  double dist = 0.0;
  for (unsigned int d = 0; d < dim; ++d)
    {
      double delta = point_2[d] - point_2[d];
      // TODO: This is poorly optimized
      for (const auto &periodic_pair : pair_set)
        {
          if (periodic_pair.direction == d)
            {
              const double length      = size[d];
              const double half_length = length / 2.0;
              delta                    = pmod(delta - half_length, length) - half_length;
            }
        }
      dist += delta * delta;
    }
  return std::sqrt(dist);
}

template <unsigned int dim>
void
RectangularMesh<dim>::declare_parameters(dealii::ParameterHandler &parameter_handler)
{
  parameter_handler.enter_subsection("Rectangular mesh");
  {
    for (const auto &dir : axis_labels)
      {
        const std::string axis {dir};

        parameter_handler.declare_entry(axis + " size",
                                        "0.0",
                                        dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                        "The size of the domain in the " + axis +
                                          "-direction.");
        parameter_handler.declare_alias(axis + " size", "size" + axis);

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
        size[i]         = parameter_handler.get_double(axis + " size");
        lower_bound[i]  = parameter_handler.get_double(axis + " lower bound");
        subdivisions[i] = parameter_handler.get_integer(axis + " subdivisions");
        if (parameter_handler.get_bool(axis + " periodic"))
          {
            periodic_directions.insert(i);
          }
      }
  }
  parameter_handler.leave_subsection();
}

template <unsigned int dim>
void
RectangularMesh<dim>::validate() const
{
  AssertThrow(size.norm() != 0.0, dealii::ExcMessage("Size of the mesh is zero."));
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
SphericalMesh<dim>::declare_parameters(dealii::ParameterHandler &parameter_handler)
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
SpatialDiscretization<dim>::declare(dealii::ParameterHandler &parameter_handler,
                                    unsigned int              n_subsections)
{
  parameter_handler.declare_entry("mesh type",
                                  "rectangular",
                                  dealii::Patterns::Selection(
                                    "rectangular|spherical|custom"),
                                  "The type of mesh to use.",
                                  true);
  SphericalMesh<dim>::declare_parameters(parameter_handler);
  RectangularMesh<dim>::declare_parameters(parameter_handler);

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

  for (unsigned int criterion_id = 0; criterion_id < n_subsections; criterion_id++)
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
                                   unsigned int              n_subsections)
{
  rectangular_mesh.assign_parameters(parameter_handler);
  spherical_mesh.assign_parameters(parameter_handler);

  global_refinement =
    (static_cast<unsigned int>(parameter_handler.get_integer("global refinement")));

  has_adaptivity = (parameter_handler.get_bool("mesh adaptivity"));

  remeshing_period =
    (static_cast<unsigned int>(parameter_handler.get_integer("remeshing period")));

  max_refinement =
    (static_cast<unsigned int>(parameter_handler.get_integer("max refinement")));
  min_refinement =
    (static_cast<unsigned int>(parameter_handler.get_integer("min refinement")));

  for (unsigned int criterion_id = 0; criterion_id < n_subsections; criterion_id++)
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
  get_mesh().validate();

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
const Mesh<dim> &
SpatialDiscretization<dim>::get_mesh() const
{
  if (mesh_type == TriangulationType::Rectangular)
    {
      return rectangular_mesh;
    }
  else if (mesh_type == TriangulationType::Spherical)
    {
      return spherical_mesh;
    }
  else if (mesh_type == TriangulationType::Custom)
    {
      AssertThrow(custom_mesh != nullptr,
                  dealii::ExcMessage("Custom mesh pointer is null."));
      return *custom_mesh;
    }
  else
    {
      Assert(false, dealii::ExcMessage("Invalid mesh type."));
      return rectangular_mesh; // This line will never be reached
    }
}

template <unsigned int dim>
Mesh<dim> &
SpatialDiscretization<dim>::get_mesh()
{
  if (mesh_type == TriangulationType::Rectangular)
    {
      return rectangular_mesh;
    }
  else if (mesh_type == TriangulationType::Spherical)
    {
      return spherical_mesh;
    }
  else if (mesh_type == TriangulationType::Custom)
    {
      AssertThrow(custom_mesh != nullptr,
                  dealii::ExcMessage("Custom mesh pointer is null."));
      return *custom_mesh;
    }
  else
    {
      Assert(false, dealii::ExcMessage("Invalid mesh type."));
      return rectangular_mesh; // This line will never be reached
    }
}

template <unsigned int dim>
void
SpatialDiscretization<dim>::generate_mesh(Triangulation &triangulation) const
{
  get_mesh().generate_mesh(triangulation);
}

template <unsigned int dim>
void
SpatialDiscretization<dim>::mark_boundaries(Triangulation &triangulation) const
{
  get_mesh().mark_boundaries(triangulation);
}

template <unsigned int dim>
std::list<PeriodicPair<dim>>
SpatialDiscretization<dim>::periodicity_set() const
{
  return get_mesh().periodicity_set();
}

template <unsigned int dim>
double
SpatialDiscretization<dim>::distance(const dealii::Point<dim> &point_1,
                                     const dealii::Point<dim> &point_2) const
{
  return get_mesh().distance(point_1, point_2);
}

template <unsigned int dim>
bool
SpatialDiscretization<dim>::should_refine_mesh(unsigned int increment) const
{
  return increment % remeshing_period == 0;
}

#include "user_inputs/spatial_discretization.inst"

PRISMS_PF_END_NAMESPACE
