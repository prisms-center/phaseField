// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef spatial_discretization_h
#define spatial_discretization_h

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/refinement_criterion.h>
#include <prismspf/utilities/utilities.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Struct that holds spatial discretization parameters.
 */
template <int dim>
struct spatialDiscretization
{
public:
  /**
   * \brief Internal enum for various triangulation types.
   */
  enum TriangulationType : std::uint8_t
  {
    rectangular,
    spherical
  };

  /**
   * \brief Constructor.
   */
  spatialDiscretization()
    : subdivisions(dim, 1) {};

  /**
   * \brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate();

  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Triangulation type
  TriangulationType type = TriangulationType::rectangular;

  // Domain extents in each cartesian direction
  dealii::Tensor<1, dim, double> size;

  // Radius of the spherical domain
  double radius = 0.0;

  // Mesh subdivisions in each cartesian direction
  std::vector<unsigned int> subdivisions;

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

template <int dim>
inline void
spatialDiscretization<dim>::postprocess_and_validate()
{
  // Assign the triangulation type
  if (radius != 0.0 && size.norm() == 0.0)
    {
      type = TriangulationType::spherical;
    }
  else if (radius == 0.0 && size.norm() != 0.0)
    {
      type = TriangulationType::rectangular;
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
          AssertThrow((criterion.value_lower_bound <= criterion.value_upper_bound),
                      dealii::ExcMessage(
                        "The lower bound of the value-based refinement "
                        "criteria must be less than or equal to the upper bound."));
        }
    }
}

template <int dim>
inline void
spatialDiscretization<dim>::print_parameter_summary() const
{
  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Spatial Discretization\n"
    << "================================================\n";

  if (type == TriangulationType::spherical)
    {
      conditionalOStreams::pout_summary() << "Domain radius: " << radius << "\n";
    }
  else if (type == TriangulationType::rectangular)
    {
      if constexpr (dim == 1)
        {
          conditionalOStreams::pout_summary()
            << "Domain size: x=" << size[0] << "\n"
            << "Subdivisions: x=" << subdivisions[0] << "\n";
        }
      else if constexpr (dim == 2)
        {
          conditionalOStreams::pout_summary()
            << "Domain size: x=" << size[0] << ", y=" << size[1] << "\n"
            << "Subdivisions: x=" << subdivisions[0] << ", y=" << subdivisions[1] << "\n";
        }
      else if constexpr (dim == 3)
        {
          conditionalOStreams::pout_summary()
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

  conditionalOStreams::pout_summary()
    << "Global refinement: " << global_refinement << "\n"
    << "Degree: " << degree << "\n"
    << "Adaptivity enabled: " << bool_to_string(has_adaptivity) << "\n"
    << "Max refinement: " << max_refinement << "\n"
    << "Min refinement: " << min_refinement << "\n"
    << "Remeshing period: " << remeshing_period << "\n";

  if (!refinement_criteria.empty())
    {
      conditionalOStreams::pout_summary() << "Refinement criteria:\n";
    }
  for (const auto &criterion : refinement_criteria)
    {
      conditionalOStreams::pout_summary()
        << "  Variable name: " << criterion.variable_name << "\n"
        << "  Variable index: " << criterion.variable_index << "\n"
        << "  Criterion type: " << criterion.criterion_to_string() << "\n"
        << "  Value lower bound: " << criterion.value_lower_bound << "\n"
        << "  Value upper bound: " << criterion.value_upper_bound << "\n"
        << "  Gradient lower bound: " << criterion.gradient_lower_bound << "\n\n";
    }
  conditionalOStreams::pout_summary() << "\n" << std::flush;
}

PRISMS_PF_END_NAMESPACE

#endif
