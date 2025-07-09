// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/refinement_criterion.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

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

private:
  // Triangulation type
  TriangulationType type = TriangulationType::Rectangular;

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