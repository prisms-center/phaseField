// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/grid_refiner_criterion.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/boundary_parameters.h>

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

template <unsigned int dim>
using Triangulation =
  std::conditional_t<dim == 1,
                     dealii::Triangulation<dim>,
                     dealii::parallel::distributed::Triangulation<dim>>;

/**
 * @brief Class for rectangular mesh parameters.
 */
template <unsigned int dim>
struct RectangularMesh
{
  /**
   * @brief Constructor.
   */
  RectangularMesh(dealii::Tensor<1, dim, double> _size,
                  std::array<unsigned int, dim>  _subdivisions)
    : size(_size)
    , subdivisions(_subdivisions) {};

  /**
   * @brief Generate the mesh.
   */
  void
  generate_mesh(Triangulation<dim> &triangulation) override
  {
    validate();
    dealii::GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                      subdivisions,
                                                      dealii::Point<dim>(),
                                                      dealii::Point<dim>(size));
    mark_boundaries(triangulation);
    mark_periodic(triangulation);
  };

  /**
   * @brief Validate
   */
  void
  validate()
  {}

  /**
   * @brief Domain extents in each cartesian direction.
   */
  dealii::Tensor<1, dim, double> size;

  /**
   * @brief Mesh subdivisions in each cartesian direction.
   */
  std::vector<unsigned int> subdivisions = std::vector<unsigned int>(dim, 1);

  /**
   * @brief Which directions have periodic conditions
   */
  std::set<unsigned int> periodic_directions;

private:
  /**
   * @brief Mark the boundaries of the mesh.
   */
  void
  mark_boundaries(Triangulation<dim> &triangulation) const
  {
    // Loop through the cells
    for (const auto &cell : triangulation->active_cell_iterators())
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
  };

  /**
   * @brief Mark the periodic faces of the mesh.
   */
  void
  mark_periodic(Triangulation<dim> &triangulation) const
  {
    // Create a vector of matched pairs that we fill and enforce upon the
    // constaints
    std::vector<
      dealii::GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
      periodicity_vector;
    for (unsigned int direction : periodic_directions)
      {
        // Grab the offset vector from one vertex to another
        dealii::Tensor<1, dim> offset;
        offset[direction] = size[direction];

        // Collect the matched pairs on the coarsest level of the mesh
        unsigned int boundary_id = direction * 2;
        dealii::GridTools::collect_periodic_faces(triangulation,
                                                  boundary_id,
                                                  boundary_id + 1,
                                                  direction,
                                                  periodicity_vector,
                                                  offset);
      }
    // Add periodicity
    triangulation.add_periodicity(periodicity_vector);
  };
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
  explicit SphericalMesh(double _radius)
    : radius(_radius)
  {}

  /**
   * @brief Generate the mesh.
   */
  void
  generate_mesh(Triangulation<dim> &triangulation) override
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
  validate()
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
public:
  /**
   * @brief Constructor.
   */
  SpatialDiscretization() = default;

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
            << "  Criterion type: " << criterion.criterion_string() << "\n"
            << "  Value lower bound: " << criterion.value_lower_bound << "\n"
            << "  Value upper bound: " << criterion.value_upper_bound << "\n"
            << "  Gradient lower bound: " << criterion.gradient_lower_bound << "\n\n";
        }
    }

  ConditionalOStreams::pout_summary() << "\n" << std::flush;
}

PRISMS_PF_END_NAMESPACE
