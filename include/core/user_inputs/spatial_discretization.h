#ifndef spatial_discretization_h
#define spatial_discretization_h

#include <core/conditional_ostreams.h>
#include <core/refinement/refinement_criterion.h>
#include <utilities/utilities.h>

/**
 * \brief Class that holds spatial discretization parameters.
 */
template <int dim>
class spatialDiscretization
{
public:
  /**
   * \brief Constructor.
   */
  spatialDiscretization()
    : subdivisions(dim, 1) {};

  /**
   * \brief Destructor.
   */
  ~spatialDiscretization() = default;

  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Domain extents in each cartesian direction
  dealii::Tensor<1, dim, double> domain_size;

  // Mesh subdivisions in each cartesian direction
  std::vector<uint> subdivisions;

  // Global refinement of mesh
  uint refine_factor = 0;

  // Element polynomial degree
  uint degree = 1;

  // Whether adaptive meshing (AMR) is enabled
  bool has_adaptivity = false;

  // Maximum global refinement for AMR
  uint max_refinement_level = 0;

  // Minimum global refinement for AMR
  uint min_refinement_level = 0;

  // The number of steps between remeshing
  uint remeshing_frequency = UINT_MAX;

  // The criteria used for remeshing
  std::vector<RefinementCriterion> refinement_criteria;
};

template <int dim>
inline void
spatialDiscretization<dim>::print_parameter_summary() const
{
  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "\tSpatial Discretization\n"
    << "================================================\n";

  if constexpr (dim == 1)
    {
      conditionalOStreams::pout_summary() << "Domain size: x=" << domain_size[0] << "\n";
    }
  else if constexpr (dim == 2)
    {
      conditionalOStreams::pout_summary()
        << "Domain size: x=" << domain_size[0] << ", y=" << domain_size[1] << "\n";
    }
  else if constexpr (dim == 3)
    {
      conditionalOStreams::pout_summary()
        << "Domain size: x=" << domain_size[0] << ", y=" << domain_size[1]
        << ", z=" << domain_size[2] << "\n";
    }

  conditionalOStreams::pout_summary()
    << "Global refinement: " << refine_factor << "\n"
    << "Degree: " << degree << "\n"
    << "Adaptivity enabled: " << bool_to_string(has_adaptivity) << "\n"
    << "Max refinement level: " << max_refinement_level << "\n"
    << "Min refinement level: " << min_refinement_level << "\n"
    << "Remeshing frequency: " << remeshing_frequency << "\n";

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

#endif