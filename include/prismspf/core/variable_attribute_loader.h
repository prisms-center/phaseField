// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <map>
#include <set>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Class to manage the variable attributes that the user specifies.
 */
class VariableAttributeLoader
{
public:
  /**
   * @brief Constructor.
   */
  VariableAttributeLoader() = default;

  /**
   * @brief Destructor.
   */
  virtual ~VariableAttributeLoader() = default;

  /**
   * @brief Copy constructor.
   */
  VariableAttributeLoader(const VariableAttributeLoader &variable_attribute_loader) =
    delete;

  /**
   * @brief Copy assignment.
   */
  VariableAttributeLoader &
  operator=(const VariableAttributeLoader &variable_attribute_loader) = delete;

  /**
   * @brief Move constructor.
   */
  VariableAttributeLoader(VariableAttributeLoader &&variable_attribute_loader) noexcept =
    delete;

  /**
   * @brief Move assignment.

   */
  VariableAttributeLoader &
  operator=(VariableAttributeLoader &&variable_attribute_loader) noexcept = delete;

  /**
   * @brief Initialize the variable attributes from the two user-facing methods
   * `load_variable_attributes()` and `loadPostProcessorVariableAttributes()`. This must
   * be called after the default constructor for derived classes.
   */
  void
  init_variable_attributes();

  // cppcheck-suppress-begin returnByReference

  /**
   * @brief getter function for variable attributes list (copy).
   *
   * This must be a copy because we don't need the class after the VariableAttributes are
   * processed.
   */
  [[nodiscard]] std::map<unsigned int, VariableAttributes>
  get_var_attributes() const;

  // cppcheck-suppress-end returnByReference

  /**
   * @brief User-implemented method where the variable attributes are set for all fields.
   */
  virtual void
  load_variable_attributes() = 0;

  /**
   * @brief Set the name of the variable at `index` to `name`.
   *
   * @param index Index of variable
   * @param name Name of variable at `index`
   */
  void
  set_variable_name(const unsigned int &index, const std::string &name);

  /**
   * @brief Set the field type of the variable at `index` to `field_type` where
   * `field_type` can be `Scalar` or `Vector`.
   *
   * @param index Index of variable
   * @param field_type Field type of variable at `index` (`Scalar` or `Vector`).
   */
  void
  set_variable_type(const unsigned int &index, const FieldInfo::TensorRank &field_type);

  /**
   * @brief Set the PDE type of the variable at `index` to `pde_type` where
   *`pde_type`can be `ExplicitTimeDependent`, `ImplicitTimeDependent`,
   *`TimeIndependent`, `Auxiliary`.
   *
   * @param index Index of variable
   * @param pde_type PDE type of variable at `index`.
   */
  void
  set_variable_equation_type(const unsigned int &index, const PDEType &pde_type);

  /**
   * @brief Set the whether the field is a postprocessed field.
   *
   * @param index Index of variable
   * @param is_postprocess Whether the field is postprocessed.
   */
  void
  set_is_postprocessed_field(const unsigned int &index, const bool &is_postprocess);

  /**
   * @brief Set whether the field is a nucleation rate.
   *
   * @param index Index of variable
   * @param is_nucleation Whether the field is a nucleation rate.
   */
  template <typename Iterable>
  void
  set_is_nucleation_rate(const unsigned int &index,
                         const bool         &is_nucleation,
                         const Iterable     &nucleating_fields);

  /**
   * @brief Set the solve block of the field.
   *
   * @param index Index of variable
   * @param solve_block Solve block.
   */
  void
  set_solve_block(const unsigned int &index, const Types::Index &solve_block);

  /**
   * @brief Add dependencies for the value term of the RHS equation of the variable at
   * `index`.
   *
   * @param index Index of variable
   * @param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_value_term_rhs(const unsigned int &index,
                                  const std::string  &dependencies);

  /**
   * @brief Add dependencies for the gradient term of the RHS equation of the variable
   * at `index`.
   *
   * @param index Index of variable
   * @param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_gradient_term_rhs(const unsigned int &index,
                                     const std::string  &dependencies);

  /**
   * @brief Add dependencies for the value term of the LHS equation of the variable at
   * `index`.
   *
   * @param index Index of variable
   * @param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_value_term_lhs(const unsigned int &index,
                                  const std::string  &dependencies);

  /**
   * @brief Add dependencies for the gradient term of the LHS equation of the variable
   * at `index`.
   *
   * @param index Index of variable
   * @param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_gradient_term_lhs(const unsigned int &index,
                                     const std::string  &dependencies);

  /**
   * @brief Insert dependencies for the value term of the RHS equation of the variable at
   * `index`.
   *
   * @param index Index of variable
   * @param dependencies Container containing list of dependency strings for
   * variable at `index` Hint: {"variable", "grad(variable)", "hess(variable)"}
   */
  template <typename Iterable>
  void
  insert_dependencies_value_term_rhs(const unsigned int &index,
                                     const Iterable     &dependencies);

  /**
   * @brief Insert dependencies for the gradient term of the RHS equation of the variable
   * at `index`.
   *
   * @param index Index of variable
   * @param dependencies Container containing list of dependency strings for
   * variable at `index` Hint: {"variable", "grad(variable)", "hess(variable)"}
   */
  template <typename Iterable>
  void
  insert_dependencies_gradient_term_rhs(const unsigned int &index,
                                        const Iterable     &dependencies);

  /**
   * @brief Insert dependencies for the value term of the LHS equation of the variable at
   * `index`.
   *
   * @param index Index of variable
   * @param dependencies Container containing list of dependency strings for
   * variable at `index` Hint: {"variable", "grad(variable)", "hess(variable)"}
   */
  template <typename Iterable>
  void
  insert_dependencies_value_term_lhs(const unsigned int &index,
                                     const Iterable     &dependencies);

  /**
   * @brief Insert dependencies for the gradient term of the LHS equation of the variable
   * at `index`.
   *
   * @param index Index of variable
   * @param dependencies Container containing list of dependency strings for
   * variable at `index` Hint: {"variable", "grad(variable)", "hess(variable)"}
   */
  template <typename Iterable>
  void
  insert_dependencies_gradient_term_lhs(const unsigned int &index,
                                        const Iterable     &dependencies);

private:
  /**
   * @brief The solutions variable & postprocessing variable attributes
   */
  std::map<unsigned int, VariableAttributes> var_attributes;

  /**
   * @brief Perform a suite of assertions on the attributes to ensure that
   * the user's inputs are well-formed
   */
  void
  validate_attributes();

  /**
   * @brief Validate that the variable name is not empty and does not contain any
   * forbidden substrings (names).
   */
  void
  validate_variable_name(const std::string           &name,
                         const std::set<std::string> &forbidden_names,
                         const std::string           &context,
                         unsigned int                 index);

  /**
   * @brief Populate dependencies that we should expect from the user.
   */
  void
  populate_dependencies(
    const std::set<std::pair<std::string, std::string>> &reg_delimiters,
    const std::set<std::pair<std::string, std::string>> &dep_type_delimiters,
    const std::string                                   &variable_name,
    unsigned int                                         index,
    std::set<std::string>                               &reg_possible_deps,
    std::map<unsigned int, std::set<std::string>>       &change_possible_deps);

  /**
   * @brief Validate the dependencies (RHS or LHS) that the user has provided.
   */
  void
  validate_dependencies(
    const std::set<std::string>                         &dependencies,
    const std::string                                   &context,
    unsigned int                                         index,
    const std::string                                   &variable_name,
    const std::set<std::string>                         &reg_possible_deps,
    const std::map<unsigned int, std::set<std::string>> &change_possible_deps);

  /**
   * @brief Validate the old solution dependencies that the user has provided.
   *
   * There are two criteria:
   * 1. The storage of old fields must be sequential (old_1 must be included before
   * old_2).
   * 2. Old fields cannot be stored for constant equations.
   */
  void
  validate_old_solution_dependencies();

  /**
   * @brief Compute the subset of VariableAttributes that belongs to a given
   * FieldSolveType and solver order.
   *
   * This function creates and returns a map of the VariablesAttributes that belong to a
   * FieldSolveType and solve order.
   *
   * @param[in] variable_attributes The set of variable attributes that we are taking the
   * subset of
   * @param[in] field_solve_type The FieldSolveType for the subset
   * @param[in] solve_priority The solve priority for the subset
   */
  [[nodiscard]] std::map<Types::Index, VariableAttributes *>
  compute_subset_attributes(
    std::map<Types::Index, VariableAttributes> &variable_attributes,
    FieldSolveType                              field_solve_type,
    Types::Index                                solve_priority) const;

  /**
   * @brief Compute the shared dependencies for a subset of VariableAttributes that belong
   * to a given FieldSolveType and solve order.
   *
   * This function computes the shared dependencies for the RHS. It updates the
   * eval_flag_set_rhs and dependency_set_rhs to be the same between all the variables.
   *
   * @param[in] variable_attributes The set of variable attributes that we are computing
   * shared dependencies for
   */
  void
  compute_shared_dependencies(
    std::map<Types::Index, VariableAttributes *> &variable_attributes);
};

PRISMS_PF_END_NAMESPACE
