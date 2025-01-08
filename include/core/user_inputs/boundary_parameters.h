#ifndef boundary_parameters_h
#define boundary_parameters_h

#include <deal.II/base/point.h>

#include <boost/algorithm/string/predicate.hpp>

#include <core/conditional_ostreams.h>
#include <core/exceptions.h>
#include <core/type_enums.h>
#include <core/variable_attributes.h>
#include <map>

/**
 * \brief Class that stores relevant information for boundary conditions of a certain
 * field.
 */
template <int dim>
class boundaryCondition
{
public:
  /**
   * \brief Constructor.
   */
  boundaryCondition() = default;

  /**
   * \brief Destructor.
   */
  ~boundaryCondition() = default;

  // Map of boundary conditions and domain boundary for which they correspond to. For a
  // simple geometry like a square the boundary ids are marked, in order, by x=0, x=max,
  // y=0, y=max. More complex geometries can have somewhat arbitrary ordering, but will
  // render some of our assertiosn moot.
  std::map<uint, boundaryType> boundary_condition_map;

  // A map of boundary values for dirichlet boundary conditions
  std::map<uint, double> dirichlet_value_map;
};

/**
 * \brief Class that holds boundary parameters.
 */
template <int dim>
class boundaryParameters
{
public:
  /**
   * \brief Constructor.
   */
  boundaryParameters() = default;

  /**
   * \brief Destructor.
   */
  ~boundaryParameters() = default;

  /**
   * \brief Compute the boundary conditions from the unfiltered string list.
   */
  void
  compute_boundary_conditions(const AttributesList &var_attributes);

  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Map of unfiltered boundary conditions strings. The first key is the global index. The
  // second key is the number of dimensions.
  std::map<uint, std::map<uint, std::string>> BC_list;

  // Map of pinned points. The first key is the global index. The pair is the pinned
  // value and point.
  std::map<uint, std::pair<double, dealii::Point<dim>>> pinned_point_list;

  // Map of boundary conditions. The first key is the global index. The second key is the
  // number of dimensions.
  std::map<uint, std::map<uint, boundaryCondition<dim>>> boundary_condition_list;
};

template <int dim>
inline void
boundaryParameters<dim>::compute_boundary_conditions(const AttributesList &var_attributes)
{
  for (const auto &[index, variable] : var_attributes)
    {
      // Ensure that boundary conditions are specified for all variables and their
      // components
      if (variable.field_type == fieldType::VECTOR)
        {
          for (uint i = 0; i < dim; i++)
            {
              AssertThrow(!BC_list.at(index).at(i).empty(),
                          dealii::ExcMessage("Boundary conditions must be specified "
                                             "for all components in all vector field."));
            }
        }
      else
        {
          AssertThrow(!BC_list.at(index).at(0).empty(),
                      dealii::ExcMessage("Boundary conditions must be specified "
                                         "for all scalar fields."));
        }
    }

  // Helper function that splits boundary condition input (delimiter ",") and assigns
  // the correct boundary condition enum.
  auto set_boundary =
    [&](const std::string &BC_string, const uint &index, const uint &component)
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
    boundaryCondition<dim> condition;
    for (uint i = 0; i < (2 * dim); i++)
      {
        if (boost::iequals(BC_string_list[i], "NATURAL"))
          {
            condition.boundary_condition_map.emplace(i, boundaryType::NATURAL);
          }
        else if (boost::iequals(BC_string_list[i].substr(0, 9), "DIRICHLET"))
          {
            condition.boundary_condition_map.emplace(i, boundaryType::DIRICHLET);
            std::string dirichlet_value =
              BC_string_list[i].substr(10, BC_string_list[i].size());
            dirichlet_value = dealii::Utilities::trim(dirichlet_value);
            condition.dirichlet_value_map.emplace(i,
                                                  dealii::Utilities::string_to_double(
                                                    dirichlet_value));
          }
        else if (boost::iequals(BC_string_list[i], "PERIODIC"))
          {
            condition.boundary_condition_map.emplace(i, boundaryType::PERIODIC);
          }
        else if (boost::iequals(BC_string_list[i].substr(0, 7), "NEUMANN"))
          {
            AssertThrow(false, FeatureNotImplemented("Neumann boundary conditions"));
          }
        else if (boost::iequals(BC_string_list[i], "NON_UNIFORM_DIRICHLET"))
          {
            condition.boundary_condition_map.emplace(i,
                                                     boundaryType::NON_UNIFORM_DIRICHLET);
          }
        else if (boost::iequals(BC_string_list[i], "NON_UNIFORM_NEUMANN"))
          {
            AssertThrow(false,
                        FeatureNotImplemented("Nonuniform neumann boundary conditions"));
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

    boundary_condition_list[component].emplace(index, condition);
  };
  // Process each boundary condition
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.field_type == fieldType::VECTOR)
        {
          for (uint i = 0; i < dim; i++)
            {
              set_boundary(BC_list.at(index).at(i), index, i);
            }
        }
      else
        {
          set_boundary(BC_list.at(index).at(0), index, 0);
        }
    }
}

template <int dim>
inline void
boundaryParameters<dim>::print_parameter_summary() const
{
  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "\tBoundary Parameters\n"
    << "================================================\n";

  conditionalOStreams::pout_summary() << "Index: ";
  for (const auto &[index, component_map] : boundary_condition_list)
    {
      for (const auto &[component, boundary_condition] : component_map)
        {
          conditionalOStreams::pout_summary() << index << "\n"
                                              << "  Component: " << component << "\n";
          for (const auto &[domain_id, boundary_type] :
               boundary_condition.boundary_condition_map)
            {
              conditionalOStreams::pout_summary()
                << "    Boundary id: " << domain_id << "    "
                << "Condition: " << to_string(boundary_type);
              if (boundary_type == boundaryType::DIRICHLET)
                {
                  conditionalOStreams::pout_summary()
                    << " = " << boundary_condition.dirichlet_value_map.at(domain_id);
                }
              conditionalOStreams::pout_summary() << "\n";
            }
        }
    }

  if (!pinned_point_list.empty())
    {
      conditionalOStreams::pout_summary() << "Pinned field index: ";
    }
  for (const auto &[index, point_value_map] : pinned_point_list)
    {
      conditionalOStreams::pout_summary()
        << index << "\n"
        << "  Value: " << point_value_map.first << "\n"
        << "  Point: " << point_value_map.second << "\n";
    }
  conditionalOStreams::pout_summary() << "\n" << std::flush;
}

#endif