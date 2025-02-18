#ifndef spatial_discretization_h
#define spatial_discretization_h

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/refinement_criterion.h>
#include <prismspf/utilities.h>

PRISMS_PF_BEGIN_NAMESPACE

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

  // Radius of the spherical domain
  double radius = 0.0;

  // Mesh subdivisions in each cartesian direction
  std::vector<unsigned int> subdivisions;

  // Global refinement of mesh
  unsigned int refine_factor = 0;

  // Element polynomial degree
  unsigned int degree = 1;

  // Whether adaptive meshing (AMR) is enabled
  bool has_adaptivity = false;

  // Maximum global refinement for AMR
  unsigned int max_refinement_level = 0;

  // Minimum global refinement for AMR
  unsigned int min_refinement_level = 0;

  // The number of steps between remeshing
  unsigned int remeshing_frequency = UINT_MAX;

  // The criteria used for remeshing
  std::vector<RefinementCriterion> refinement_criteria;
};

template <int dim>
inline void
spatialDiscretization<dim>::print_parameter_summary() const
{
  prisms::conditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Spatial Discretization\n"
    << "================================================\n";

  if (radius != 0.0)
    {
      prisms::conditionalOStreams::pout_summary() << "Domain radius: " << radius << "\n";
    }
  else
    {
      if constexpr (dim == 1)
        {
          prisms::conditionalOStreams::pout_summary()
            << "Domain size: x=" << domain_size[0] << "\n"
            << "Subdivisions: x=" << subdivisions[0] << "\n";
        }
      else if constexpr (dim == 2)
        {
          prisms::conditionalOStreams::pout_summary()
            << "Domain size: x=" << domain_size[0] << ", y=" << domain_size[1] << "\n"
            << "Subdivisions: x=" << subdivisions[0] << ", y=" << subdivisions[1] << "\n";
        }
      else if constexpr (dim == 3)
        {
          prisms::conditionalOStreams::pout_summary()
            << "Domain size: x=" << domain_size[0] << ", y=" << domain_size[1]
            << ", z=" << domain_size[2] << "\n"
            << "Subdivisions: x=" << subdivisions[0] << ", y=" << subdivisions[1]
            << ", z=" << subdivisions[2] << "\n";
        }
    }

  prisms::conditionalOStreams::pout_summary()
    << "Global refinement: " << refine_factor << "\n"
    << "Degree: " << degree << "\n"
    << "Adaptivity enabled: " << bool_to_string(has_adaptivity) << "\n"
    << "Max refinement level: " << max_refinement_level << "\n"
    << "Min refinement level: " << min_refinement_level << "\n"
    << "Remeshing frequency: " << remeshing_frequency << "\n";

  if (!refinement_criteria.empty())
    {
      prisms::conditionalOStreams::pout_summary() << "Refinement criteria:\n";
    }
  for (const auto &criterion : refinement_criteria)
    {
      prisms::conditionalOStreams::pout_summary()
        << "  Variable name: " << criterion.variable_name << "\n"
        << "  Variable index: " << criterion.variable_index << "\n"
        << "  Criterion type: " << criterion.criterion_to_string() << "\n"
        << "  Value lower bound: " << criterion.value_lower_bound << "\n"
        << "  Value upper bound: " << criterion.value_upper_bound << "\n"
        << "  Gradient lower bound: " << criterion.gradient_lower_bound << "\n\n";
    }
  prisms::conditionalOStreams::pout_summary() << "\n" << std::flush;
}

PRISMS_PF_END_NAMESPACE

#endif