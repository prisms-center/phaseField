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
 * @brief Base class for mesh parameters with a generator functions.
 */
template <unsigned int dim>
class Mesh
{
public:
  using Triangulation =
    std::conditional_t<dim == 1,
                       dealii::Triangulation<dim>,
                       dealii::parallel::distributed::Triangulation<dim>>;

  /**
   * @brief Constructor.
   */
  Mesh() = default;

  /**
   * @brief Virtual destructor.
   */
  virtual ~Mesh() = default;

  /**
   * @brief Generate the mesh.
   */
  virtual void
  generate_mesh([[maybe_unused]] std::shared_ptr<Triangulation> triangulation)
  {
    AssertThrow(
      is_initialized,
      dealii::ExcMessage(
        "Mesh parameters not initialized correctly. You tried to generate a mesh that "
        "has default parameters. Typically this yields a mesh with size 0.0."));
    AssertThrow(
      triangulation->n_cells() == 0,
      dealii::ExcMessage(
        "Mesh triangulation is not empty. This is likely because you "
        "tried to generate a mesh without setting the parameters first. "
        "Please ensure that you have set the parameters before generating "
        "the mesh. Another possibility is that you tried to generate multiple meshes. "
        "Check that you only have one subsection in the parameters file for the mesh."));
    is_generated = true;
  };

  /**
   * @brief Mark the boundaries of the mesh.
   */
  void
  mark_boundaries(
    [[maybe_unused]] std::shared_ptr<Triangulation> triangulation,
    [[maybe_unused]] const typename BoundaryParameters<dim>::BoundaryConditionMap
      &boundary_condition_list)
  {
    AssertThrow(is_generated,
                dealii::ExcMessage(
                  "Mesh must be generated before marking boundary ids."));
    boundaries_marked = true;
  };

  /**
   * @brief Mark the periodic faces of the mesh.
   */
  void
  mark_periodic(
    [[maybe_unused]] std::shared_ptr<Triangulation> triangulation,
    [[maybe_unused]] const typename BoundaryParameters<dim>::BoundaryConditionMap
      &boundary_condition_list)
  {
    AssertThrow(is_generated,
                dealii::ExcMessage(
                  "Mesh must be generated before marking periodic boundaries."));
    AssertThrow(boundaries_marked,
                dealii::ExcMessage(
                  "Mesh boundaries must be marked before marking periodic boundaries."));
  };

  /**
   * @brief Set that the mesh parameters are initialized correctly.
   */
  void
  set_initialized()
  {
    is_initialized = true;
  };

private:
  /**
   * @brief Whether the class is initialized correctly.
   *
   * This is necessary because the parameters are determined at runtime, leading to the
   * generation of a bunch of empty mesh objects.
   */
  bool is_initialized = false;

  /**
   * @brief Whether the mesh has been generated.
   */
  bool is_generated = false;

  /**
   * @brief Whether the boundaries have been marked.
   */
  bool boundaries_marked = false;
};

/**
 * @brief Class for rectangular mesh parameters.
 */
template <unsigned int dim>
class RectangularMesh : public Mesh<dim>
{
  /**
   * @brief Constructor.
   */
  RectangularMesh(dealii::Tensor<1, dim, double> _size,
                  std::vector<unsigned int>      _subdivisions)
    : size(_size)
  {
    Assert(_subdivisions.size() == dim,
           dealii::ExcMessage(
             "Subdivisions vector size must match the number of dimensions"));
    subdivisions = _subdivisions;

    // If the x direction is greater than 0, we set the mesh as initialized.
    // TODO (landinjm): Check that the other directions are also greater than 0.
    if (size[0] > 0.0)
      {
        this->Mesh<dim>::set_initialized();
      }
  };

  /**
   * @brief Generate the mesh.
   */
  void
  generate_mesh(std::shared_ptr<typename Mesh<dim>::Triangulation> triangulation) override
  {
    // Call base class generate mesh to check initialization
    this->Mesh<dim>::generate_mesh(triangulation);

    dealii::GridGenerator::subdivided_hyper_rectangle(*triangulation,
                                                      subdivisions,
                                                      dealii::Point<dim>(),
                                                      dealii::Point<dim>(size));
  };

  /**
   * @brief Mark the boundaries of the mesh.
   */
  void
  mark_boundaries(std::shared_ptr<typename Mesh<dim>::Triangulation> triangulation,
                  const typename BoundaryParameters<dim>::BoundaryConditionMap
                    &boundary_condition_list) override
  {
    // Call base class mark boundaries to check that the mesh has been generated
    this->Mesh<dim>::mark_boundaries(triangulation, boundary_condition_list);

    // Check that the user has not specified extra boundary conditions
    for (const auto &[index, boundary_conditions] : boundary_condition_list)
      {
        AssertThrow(boundary_conditions.get_boundary_condition_map().size() == (2 * dim),
                    dealii::ExcMessage(
                      "Rectangular meshes only have 2*dim boundaries. Please ensure that "
                      "you have not specified extra boundary conditions."));
      }

    // Loop through the cells
    for (const auto &cell : triangulation->active_cell_iterators())
      {
        // Mark the faces (faces_per_cell = 2*dim)
        for (unsigned int face_number = 0;
             face_number < dealii::GeometryInfo<dim>::faces_per_cell;
             ++face_number)
          {
            // Direction for quad and hex cells
            auto direction = static_cast<unsigned int>(std::floor(face_number / 2));

            // Mark the boundary id for x=0, y=0, z=0 and x=max, y=max, z=max
            if (std::fabs(cell->face(face_number)->center()(direction) - 0) <
                  Defaults::mesh_tolerance ||
                std::fabs(cell->face(face_number)->center()(direction) -
                          (size[direction])) < Defaults::mesh_tolerance)
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
  mark_periodic(std::shared_ptr<typename Mesh<dim>::Triangulation> triangulation,
                const typename BoundaryParameters<dim>::BoundaryConditionMap
                  &boundary_condition_list) override
  {
    // Call base class mark boundaries to check that the mesh has been generated and
    // boundaries marked
    this->Mesh<dim>::mark_periodic(triangulation, boundary_condition_list);

    // TODO (landinjm): Add some assertions here

    // Add periodicity in the triangulation where specified in the boundary conditions.
    // Note
    // that if one field is periodic all others should be as well.
    for (const auto &[index, boundary_conditions] : boundary_condition_list)
      {
        for (const auto &[component, condition] : boundary_conditions)
          {
            for (const auto &[boundary_id, boundary_type] :
                 condition.get_boundary_condition_map())
              {
                if (boundary_type == BoundaryCondition::Type::Periodic)
                  {
                    // Skip boundary ids that are odd since those map to the even faces
                    if (boundary_id % 2 != 0)
                      {
                        continue;
                      }

                    // Create a vector of matched pairs that we fill and enforce upon the
                    // constaints
                    std::vector<dealii::GridTools::PeriodicFacePair<
                      typename Mesh<dim>::Triangulation::cell_iterator>>
                      periodicity_vector;

                    // Determine the direction
                    const auto direction =
                      static_cast<unsigned int>(std::floor(boundary_id / dim));

                    // Collect the matched pairs on the coarsest level of the mesh
                    dealii::GridTools::collect_periodic_faces(*triangulation,
                                                              boundary_id,
                                                              boundary_id + 1,
                                                              direction,
                                                              periodicity_vector);

                    // Set constraints
                    triangulation->add_periodicity(periodicity_vector);
                  }
              }
          }
      }
  };

  /**
   * @brief Get the size in a certain direction.
   */
  [[nodiscard]] double
  get_size(unsigned int direction) const
  {
    Assert(direction < dim,
           dealii::ExcMessage("Direction must be less than the number of dimensions"));
    Assert(size[direction] > 0.0,
           dealii::ExcMessage("Size in the specified direction must be greater than zero "
                              "for rectangular meshes"));
    return size[direction];
  };

  /**
   * @brief Get the subdivisions in a certain direction.
   */
  [[nodiscard]] unsigned int
  get_subdivisions(unsigned int direction) const
  {
    Assert(direction < dim,
           dealii::ExcMessage("Direction must be less than the number of dimensions"));
    Assert(subdivisions[direction] > 0,
           dealii::ExcMessage(
             "Subdivisions in the specified direction must be greater than zero "
             "for rectangular meshes"));
    return subdivisions[direction];
  }

private:
  /**
   * @brief Domain extents in each cartesian direction.
   */
  dealii::Tensor<1, dim, double> size;

  /**
   * @brief Mesh subdivisions in each cartesian direction.
   */
  std::vector<unsigned int> subdivisions = std::vector<unsigned int>(dim, 1);
};

/**
 * @brief Class for spherical mesh parameters.
 */
template <unsigned int dim>
class SphericalMesh
{
public:
  /**
   * @brief Constructor.
   */
  explicit SphericalMesh(double _radius)
    : radius(_radius)
  {
    if (radius > 0.0)
      {
        this->Mesh<dim>::set_initialized();
      }
  };

  /**
   * @brief Generate the mesh.
   */
  void
  generate_mesh(std::shared_ptr<typename Mesh<dim>::Triangulation> triangulation) override
  {
    // Call base class generate mesh to check initialization
    this->Mesh<dim>::generate_mesh(triangulation);

    AssertThrow(dim != 1, dealii::ExcMessage("Spherical mesh not valid in 1D"));
    dealii::GridGenerator::hyper_ball(*triangulation, dealii::Point<dim>(), radius);
  };

  /**
   * @brief Mark the boundaries of the mesh.
   */
  void
  mark_boundaries(std::shared_ptr<typename Mesh<dim>::Triangulation> triangulation,
                  const typename BoundaryParameters<dim>::BoundaryConditionMap
                    &boundary_condition_list) override
  {
    // Call base class mark boundaries to check that the mesh has been generated
    this->Mesh<dim>::mark_boundaries(triangulation, boundary_condition_list);

    for (const auto &[index, boundary_conditions] : boundary_condition_list)
      {
        // There are no special boundary ids to mark for a spherical mesh. Check that the
        // user has specified extra boundary conditions.
        AssertThrow(boundary_conditions.get_boundary_condition_map().size() == 1,
                    dealii::ExcMessage("Spherical meshes only have one boundary. Please "
                                       "ensure that you have not specified "
                                       "extra boundary conditions."));
      }
  };

  /**
   * @brief Mark the periodic faces of the mesh.
   */
  void
  mark_periodic(std::shared_ptr<typename Mesh<dim>::Triangulation> triangulation,
                const typename BoundaryParameters<dim>::BoundaryConditionMap
                  &boundary_condition_list) override
  {
    // Call base class mark boundaries to check that the mesh has been generated and
    // boundaries marked
    this->Mesh<dim>::mark_periodic(triangulation, boundary_condition_list);

    // There are no periodic boundaries for a spherical mesh. Check that the user has not
    // specified any.
    for (const auto &[index, boundary_conditions] : boundary_condition_list)
      {
        AssertThrow(boundary_conditions.get_boundary_condition_map().begin()->second !=
                      BoundaryCondition::Type::Periodic,
                    dealii::ExcMessage("Spherical meshes cannot have periodic boundary "
                                       "conditions."));
      }
  };

  /**
   * @brief Get the radius of the spherical domain.
   */
  [[nodiscard]] double
  get_radius() const
  {
    Assert(radius > 0.0,
           dealii::ExcMessage("Spherical mesh radius must be greater than zero."));
    return radius;
  }

private:
  /**
   * @brief Radius of the spherical domain.
   */
  double radius = 0.0;
};

/**
 * @brief Struct that holds spatial discretization parameters.
 */
template <unsigned int dim>
struct SpatialDiscretization
{
public:
  /**
   * @brief Internal enum for various triangulation types.
   */
  enum TriangulationType : std::uint8_t
  {
    Rectangular,
    Spherical
  };

  /**
   * @brief Constructor.
   */
  SpatialDiscretization() = default;

  /**
   * @brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate();

  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  /**
   * @brief Whether the provided increment is a valid grid refinement step.
   */
  [[nodiscard]] bool
  should_refine_mesh(unsigned int increment) const
  {
    return increment % remeshing_period == 0;
  }

  /**
   * @brief Get the domain extents in each cartesian direction
   */
  [[nodiscard]] const dealii::Tensor<1, dim, double> &
  get_size() const
  {
    return size;
  }

  /**
   * @brief Set the domain extents in each cartesian direction
   */
  void
  set_size(const unsigned int &direction, const double &_size)
  {
    size[direction] = _size;
  }

  /**
   * @brief Get the radius of the Spherical domain
   */
  [[nodiscard]] double
  get_radius() const
  {
    return radius;
  }

  /**
   * @brief Set the radius of the Spherical domain
   */
  void
  set_radius(const double &_radius)
  {
    radius = _radius;
  }

  /**
   * @brief Get the mesh subdivisions in each cartesian direction
   */
  [[nodiscard]] const std::vector<unsigned int> &
  get_subdivisions() const
  {
    return subdivisions;
  }

  /**
   * @brief Set the mesh subdivisions in each cartesian direction
   */
  void
  set_subdivisions(const unsigned int &direction, const unsigned int &_subdivisions)
  {
    subdivisions[direction] = _subdivisions;
  }

  /**
   * @brief Get the global refinement of mesh
   */
  [[nodiscard]] unsigned int
  get_global_refinement() const
  {
    return global_refinement;
  }

  /**
   * @brief Set the global refinement of mesh
   */
  void
  set_global_refinement(const unsigned int &_global_refinement)
  {
    global_refinement = _global_refinement;
  }

  /**
   * @brief Get the element polynomial degree
   */
  [[nodiscard]] unsigned int
  get_degree() const
  {
    return degree;
  }

  /**
   * @brief Set the element polynomial degree
   */
  void
  set_degree(const unsigned int &_degree)
  {
    degree = _degree;
  }

  /**
   * @brief Get whether adaptive meshing (AMR) is enabled
   */
  [[nodiscard]] bool
  get_has_adaptivity() const
  {
    return has_adaptivity;
  }

  /**
   * @brief Set whether adaptive meshing (AMR) is enabled
   */
  void
  set_has_adaptivity(const bool &_has_adaptivity)
  {
    has_adaptivity = _has_adaptivity;
  }

  /**
   * @brief Get the maximum global refinement for AMR
   */
  [[nodiscard]] unsigned int
  get_max_refinement() const
  {
    return max_refinement;
  }

  /**
   * @brief Set the maximum global refinement for AMR
   */
  void
  set_max_refinement(const unsigned int &_max_refinement)
  {
    max_refinement = _max_refinement;
  }

  /**
   * @brief Get the minimum global refinement for AMR
   */
  [[nodiscard]] unsigned int
  get_min_refinement() const
  {
    return min_refinement;
  }

  /**
   * @brief Set the minimum global refinement for AMR
   */
  void
  set_min_refinement(const unsigned int &_min_refinement)
  {
    min_refinement = _min_refinement;
  }

  /**
   * @brief Get the number of steps between remeshing
   */
  [[nodiscard]] unsigned int
  get_remeshing_period() const
  {
    return remeshing_period;
  }

  /**
   * @brief Set the number of steps between remeshing
   */
  void
  set_remeshing_period(const unsigned int &_remeshing_period)
  {
    remeshing_period = _remeshing_period;
  }

  /**
   * @brief Get the refinement criteria
   */
  [[nodiscard]] const std::vector<GridRefinement::RefinementCriterion> &
  get_refinement_criteria() const
  {
    return refinement_criteria;
  }

  /**
   * @brief Set the refinement criteria
   */
  void
  add_refinement_criteria(
    const GridRefinement::RefinementCriterion &_refinement_criterion)
  {
    refinement_criteria.push_back(_refinement_criterion);
  }

  /**
   * @brief Get the triangulation type
   */
  [[nodiscard]] TriangulationType
  get_type() const
  {
    return type;
  }

  /**
   * @brief Set the lower bound in each cartesian direction
   */
  void
  set_lower_bound(const unsigned int &direction, const double &_lower_bound)
  {
    lower_bound[direction] = _lower_bound;
  }

  /**
   * @brief Get the lower bound for the rectangular grid
   */
  [[nodiscard]] const dealii::Point<dim, double> &
  get_lower_bound() const
  {
    return lower_bound;
  }

  /**
   * @brief Get the upper bound for the rectangular grid
   */
  [[nodiscard]] const dealii::Point<dim, double> &
  get_upper_bound() const
  {
    upper_bound = dealii::Point<dim>(size) + lower_bound;
    return upper_bound;
  }

private:
  // Triangulation type
  TriangulationType type = TriangulationType::Rectangular;

  // Domain lower bound for rectangular domains
  dealii::Point<dim, double> lower_bound;

  // Domain upper bound for rectangular domains
  mutable dealii::Point<dim, double> upper_bound;

  // Domain extents in each cartesian direction
  dealii::Tensor<1, dim, double> size;

  // Radius of the Spherical domain
  double radius = 0.0;

  // Mesh subdivisions in each cartesian direction
  std::vector<unsigned int> subdivisions = std::vector<unsigned int>(dim, 1);

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
  std::vector<GridRefinement::RefinementCriterion> refinement_criteria;
};

template <unsigned int dim>
inline void
SpatialDiscretization<dim>::postprocess_and_validate()
{
  // Assign the triangulation type
  if (radius != 0.0 && size.norm() == 0.0)
    {
      type = TriangulationType::Spherical;
    }
  else if (radius == 0.0 && size.norm() != 0.0)
    {
      type = TriangulationType::Rectangular;
    }
  else
    {
      AssertThrow(false, UnreachableCode());
    }

  // Check that AMR is not enabled for 1D
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

      // Check that the refinement criteria are valid for the lower and upper bounds
      for (const auto &criterion : refinement_criteria)
        {
          AssertThrow((criterion.get_value_lower_bound() <=
                       criterion.get_value_upper_bound()),
                      dealii::ExcMessage(
                        "The lower bound of the value-based refinement "
                        "criteria must be less than or equal to the upper bound."));
        }
    }
}

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
      ConditionalOStreams::pout_summary() << "Domain radius: " << radius << "\n";
    }
  else if (type == TriangulationType::Rectangular)
    {
      if constexpr (dim == 1)
        {
          ConditionalOStreams::pout_summary()
            << "Domain size: x=" << size[0] << "\n"
            << "Subdivisions: x=" << subdivisions[0] << "\n";
        }
      else if constexpr (dim == 2)
        {
          ConditionalOStreams::pout_summary()
            << "Domain size: x=" << size[0] << ", y=" << size[1] << "\n"
            << "Subdivisions: x=" << subdivisions[0] << ", y=" << subdivisions[1] << "\n";
        }
      else if constexpr (dim == 3)
        {
          ConditionalOStreams::pout_summary()
            << "Domain size: x=" << size[0] << ", y=" << size[1] << ", z=" << size[2]
            << "\n"
            << "Subdivisions: x=" << subdivisions[0] << ", y=" << subdivisions[1]
            << ", z=" << subdivisions[2] << "\n";
        }
    }
  else
    {
      AssertThrow(false, UnreachableCode());
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
    }
  for (const auto &criterion : refinement_criteria)
    {
      ConditionalOStreams::pout_summary()
        << "  Criterion type: " << criterion.criterion_to_string() << "\n"
        << "  Value lower bound: " << criterion.get_value_lower_bound() << "\n"
        << "  Value upper bound: " << criterion.get_value_upper_bound() << "\n"
        << "  Gradient lower bound: " << criterion.get_gradient_lower_bound() << "\n\n";
    }
  ConditionalOStreams::pout_summary() << "\n" << std::flush;
}

PRISMS_PF_END_NAMESPACE
