#ifndef boundary_parameters_h
#define boundary_parameters_h

#include <deal.II/base/point.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <boost/algorithm/string/predicate.hpp>

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/types.h>

#include <cstddef>
#include <map>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Class that stores relevant information for boundary conditions of a certain
 * field.
 */
struct boundaryCondition
{
public:
  /**
   * \brief Type of boundary condition
   */
  enum type : std::uint8_t
  {
    UNDEFINED_BOUNDARY,
    NATURAL,
    DIRICHLET,
    PERIODIC,
    NEUMANN,
    NON_UNIFORM_DIRICHLET,
    NON_UNIFORM_NEUMANN
  };

  // Map of boundary conditions and domain boundary for which they correspond to. For a
  // simple geometry like a square the boundary ids are marked, in order, by x=0, x=max,
  // y=0, y=max. More complex geometries can have somewhat arbitrary ordering, but will
  // render some of our assertions moot.
  std::map<dealii::types::boundary_id, type> boundary_condition_map;

  // A map of boundary values for dirichlet boundary conditions
  std::map<dealii::types::boundary_id, double> dirichlet_value_map;

  /**
   * \brief Enum to string for type
   */
  [[nodiscard]] std::string
  to_string(type i) const
  {
    switch (i)
      {
        case type::UNDEFINED_BOUNDARY:
          return "UNDEFINED_BOUNDARY";
        case type::NATURAL:
          return "NATURAL";
        case type::DIRICHLET:
          return "DIRICHLET";
        case type::PERIODIC:
          return "PERIODIC";
        case type::NEUMANN:
          return "NEUMANN";
        case type::NON_UNIFORM_DIRICHLET:
          return "NON_UNIFORM_DIRICHLET";
        case type::NON_UNIFORM_NEUMANN:
          return "NON_UNIFORM_NEUMANN";
        default:
          return "UNKNOWN";
      }
  }
};

/**
 * \brief Class that holds boundary parameters.
 */
template <int dim>
struct boundaryParameters
{
public:
  using BoundaryConditionMap =
    std::map<types::index, std::map<unsigned int, boundaryCondition>>;
  using BCList         = std::map<types::index, std::map<unsigned int, std::string>>;
  using PinnedPointMap = std::map<types::index, std::pair<double, dealii::Point<dim>>>;

  /**
   * \brief Compute the boundary conditions from the unfiltered string list.
   */
  void
  compute_boundary_conditions(
    const std::map<unsigned int, variableAttributes> &var_attributes);

  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Map of unfiltered boundary conditions strings. The first key is the global index. The
  // second key is the number of dimensions.
  BCList BC_list;

  // Map of pinned points. The first key is the global index. The pair is the pinned
  // value and point.
  PinnedPointMap pinned_point_list;

  // Map of boundary conditions. The first key is the global index. The second key is the
  // number of dimensions.
  BoundaryConditionMap boundary_condition_list;

private:
  /**
   * \brief Set the boundary for a single component of a field index.
   */
  void
  set_boundary(const std::string  &BC_string,
               const types::index &index,
               const unsigned int &component);

  /**
   * \brief Perform a check on the boundary conditions to ensure that they are valid
   */
  void
  validate_boundary_conditions() const;
};

template <int dim>
inline void
boundaryParameters<dim>::compute_boundary_conditions(
  const std::map<unsigned int, variableAttributes> &var_attributes)
{
  for (const auto &[index, variable] : var_attributes)
    {
      // Ensure that boundary conditions are specified for all variables and their
      // components
      if (variable.field_type == fieldType::VECTOR)
        {
          for (unsigned int i = 0; i < dim; i++)
            {
              AssertThrow(!BC_list.at(variable.field_index).at(i).empty(),
                          dealii::ExcMessage("Boundary conditions must be specified "
                                             "for all components in all vector field."));

              set_boundary(BC_list.at(variable.field_index).at(i),
                           variable.field_index,
                           i);
            }
        }
      else
        {
          AssertThrow(!BC_list.at(variable.field_index).at(0).empty(),
                      dealii::ExcMessage("Boundary conditions must be specified "
                                         "for all scalar fields."));

          set_boundary(BC_list.at(variable.field_index).at(0), variable.field_index, 0);
        }
    }

  // Validate boundary conditions
  validate_boundary_conditions();

  // Clear the BC_list now that it's no longer necessary
  BC_list.clear();
}

template <int dim>
inline void
boundaryParameters<dim>::print_parameter_summary() const
{
  prisms::conditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Boundary Parameters\n"
    << "================================================\n";

  for (const auto &[index, component_map] : boundary_condition_list)
    {
      prisms::conditionalOStreams::pout_summary() << "Index: " << index << "\n";
      for (const auto &[component, boundary_condition] : component_map)
        {
          prisms::conditionalOStreams::pout_summary()
            << "  Component: " << component << "\n";
          for (const auto &[domain_id, boundary_type] :
               boundary_condition.boundary_condition_map)
            {
              prisms::conditionalOStreams::pout_summary()
                << "    Boundary id: " << domain_id << "    "
                << "Condition: " << boundary_condition.to_string(boundary_type);
              if (boundary_type == boundaryCondition::type::DIRICHLET)
                {
                  prisms::conditionalOStreams::pout_summary()
                    << " = " << boundary_condition.dirichlet_value_map.at(domain_id);
                }
              prisms::conditionalOStreams::pout_summary() << "\n";
            }
        }
    }

  if (!pinned_point_list.empty())
    {
      prisms::conditionalOStreams::pout_summary() << "Pinned field index: ";
    }
  for (const auto &[index, point_value_map] : pinned_point_list)
    {
      prisms::conditionalOStreams::pout_summary()
        << index << "\n"
        << "  Value: " << point_value_map.first << "\n"
        << "  Point: " << point_value_map.second << "\n";
    }
  prisms::conditionalOStreams::pout_summary() << "\n" << std::flush;
}

template <int dim>
inline void
boundaryParameters<dim>::set_boundary(const std::string  &BC_string,
                                      const types::index &index,
                                      const unsigned int &component)
{
  // Split string
  auto BC_string_list = dealii::Utilities::split_string_list(BC_string);

  // Check that there is either 1 or 2*dim entries in the vector. This can be changed
  // later to support other geometries.
  AssertThrow(
    BC_string_list.size() == 1 || BC_string_list.size() == static_cast<size_t>(2 * dim),
    dealii::ExcMessage("Either 1 or 2*dim boundary conditions must be specified."));

  // If there is only 1 boundary condition resize BC_string_list, copying the first
  // entry.
  if (BC_string_list.size() == 1)
    {
      BC_string_list.resize(static_cast<size_t>(2 * dim), BC_string_list[0]);
    }

  // Assign boundary condition
  boundaryCondition condition;
  for (unsigned int i = 0; i < (2 * dim); i++)
    {
      if (boost::iequals(BC_string_list[i], "NATURAL"))
        {
          condition.boundary_condition_map.emplace(i, boundaryCondition::type::NATURAL);
        }
      else if (boost::iequals(BC_string_list[i].substr(0, 9), "DIRICHLET"))
        {
          condition.boundary_condition_map.emplace(i, boundaryCondition::type::DIRICHLET);
          std::string dirichlet_value =
            BC_string_list[i].substr(10, BC_string_list[i].size());
          dirichlet_value = dealii::Utilities::trim(dirichlet_value);
          condition.dirichlet_value_map.emplace(i,
                                                dealii::Utilities::string_to_double(
                                                  dirichlet_value));
        }
      else if (boost::iequals(BC_string_list[i], "PERIODIC"))
        {
          condition.boundary_condition_map.emplace(i, boundaryCondition::type::PERIODIC);
        }
      else if (boost::iequals(BC_string_list[i].substr(0, 7), "NEUMANN"))
        {
          AssertThrow(false,
                      prisms::FeatureNotImplemented("Neumann boundary conditions"));
        }
      else if (boost::iequals(BC_string_list[i], "NON_UNIFORM_DIRICHLET"))
        {
          condition.boundary_condition_map
            .emplace(i, boundaryCondition::type::NON_UNIFORM_DIRICHLET);
        }
      else if (boost::iequals(BC_string_list[i], "NON_UNIFORM_NEUMANN"))
        {
          AssertThrow(false,
                      prisms::FeatureNotImplemented(
                        "Nonuniform neumann boundary conditions"));
        }
      else
        {
          AssertThrow(false,
                      dealii::ExcMessage("Invalid boundary condition " +
                                         BC_string_list[i]));
        }
      // If periodic boundary conditions are used, ensure that they are applied on
      // both sides of the domain.
      if (i % 2 == 0)
        {
          AssertThrow(boost::iequals(BC_string_list[i], "PERIODIC") ==
                        boost::iequals(BC_string_list[i + 1], "PERIODIC"),
                      dealii::ExcMessage("Periodic boundary condition must be "
                                         "specified on both sides of domain"));
        }
    }

  boundary_condition_list[index].emplace(component, condition);
}

template <int dim>
inline void
boundaryParameters<dim>::validate_boundary_conditions() const
{
  // Throw a warning if the pinned point is not on a vertex
  for (const auto &[index, point_value_map] : pinned_point_list)
    {
      const auto point     = point_value_map.second;
      bool       on_vertex = false;
      if constexpr (dim == 1)
        {
          const dealii::Point<1> vertex_1(0);

          on_vertex = point == vertex_1;
        }
      else if constexpr (dim == 2)
        {
          const dealii::Point<2> vertex_1(0, 0);

          on_vertex = point == vertex_1;
        }
      else if constexpr (dim == 3)
        {
          const dealii::Point<3> vertex_1(0, 0, 0);

          on_vertex = point == vertex_1;
        }
      else
        {
          AssertThrow(false, prisms::UnreachableCode());
        }

      AssertThrow(on_vertex, dealii::ExcMessage("Pinned point must be on the origin"));
    }

  // Throw a warning if only some fields have periodic boundary conditions
  std::vector<bool> periodic_ids(static_cast<dealii::types::boundary_id>(2 * dim), false);
  for (const auto &[index, component_map] : boundary_condition_list)
    {
      for (const auto &[component, boundary_condition] : component_map)
        {
          for (const auto &[domain_id, boundary_type] :
               boundary_condition.boundary_condition_map)
            {
              if (boundary_type == boundaryCondition::type::PERIODIC)
                {
                  periodic_ids[domain_id] = true;
                }
            }
        }
    }
  for (const auto &[index, component_map] : boundary_condition_list)
    {
      for (const auto &[component, boundary_condition] : component_map)
        {
          for (const auto &[domain_id, boundary_type] :
               boundary_condition.boundary_condition_map)
            {
              if (boundary_type != boundaryCondition::type::PERIODIC &&
                  periodic_ids[domain_id])
                {
                  AssertThrow(
                    false,
                    dealii::ExcMessage(
                      "All fields for a given domain id (side) must have periodic "
                      "boundary conditions if any field has periodic boundary "
                      "conditions"));
                }
            }
        }
    }
}

PRISMS_PF_END_NAMESPACE

#endif