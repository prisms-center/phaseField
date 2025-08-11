#ifndef variable_attribute_loader_h
#define variable_attribute_loader_h

#include <prismspf/config.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <map>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Class to manage the variable attributes that the user specifies.
 */
class variableAttributeLoader
{
public:
  /**
   * \brief Constructor.
   */
  variableAttributeLoader() = default;

  /**
   * \brief Destructor.
   */
  virtual ~variableAttributeLoader() = default;

  /**
   * \brief Initialize the variable attributes from the two user-facing methods
   * `loadVariableAttributes()` and `loadPostProcessorVariableAttributes()`. This must be
   * called after the default constructor for derived classes.
   */
  void
  init_variable_attributes();

  /**
   * \brief getter function for variable attributes list (copy).
   */
  [[nodiscard]] std::map<unsigned int, variableAttributes>
  get_var_attributes() const;

protected:
  /**
   * \brief User-facing method where the variable attributes are set for solution fields
   * and postprocess fields.
   */
  virtual void
  loadVariableAttributes();

  /**
   * \brief Set the name of the variable at `index` to `name`.
   *
   * \param index Index of variable
   * \param name Name of variable at `index`
   */
  void
  set_variable_name(const unsigned int &index, const std::string &name);

  /**
   * \brief Set the field type of the variable at `index` to `field_type` where
   * `field_type` can be `SCALAR` or `VECTOR`.
   *
   * \param index Index of variable
   * \param field_type Field type of variable at `index` (`SCALAR` or `VECTOR`).
   */
  void
  set_variable_type(const unsigned int &index, const fieldType &field_type);

  /**
   * \brief Set the PDE type of the variable at `index` to `pde_type` where
   *`pde_type`can be `EXPLICIT_TIME_DEPENDENT`, `IMPLICIT_TIME_DEPENDENT`,
   *`TIME_INDEPENDENT`, `AUXILIARY`.
   *
   * \param index Index of variable
   * \param pde_type PDE type of variable at `index`.
   */
  void
  set_variable_equation_type(const unsigned int &index, const PDEType &pde_type);

  /**
   * \brief Set the whether the field is a postprocessed field.
   *
   * \param index Index of variable
   * \param is_postprocess Whether the field is postprocessed.
   */
  void
  set_is_postprocessed_field(const unsigned int &index, const bool &is_postprocess);

  /**
   * \brief Add dependencies for the value term of the RHS equation of the variable at
   * `index`.
   *
   * \param index Index of variable
   * \param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_value_term_RHS(const unsigned int &index,
                                  const std::string  &dependencies);

  /**
   * \brief Add dependencies for the gradient term of the RHS equation of the variable
   * at `index`.
   *
   * \param index Index of variable
   * \param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_gradient_term_RHS(const unsigned int &index,
                                     const std::string  &dependencies);

  /**
   * \brief Add dependencies for the value term of the LHS equation of the variable at
   * `index`.
   *
   * \param index Index of variable
   * \param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_value_term_LHS(const unsigned int &index,
                                  const std::string  &dependencies);

  /**
   * \brief Add dependencies for the gradient term of the LHS equation of the variable
   * at `index`.
   *
   * \param index Index of variable
   * \param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_gradient_term_LHS(const unsigned int &index,
                                     const std::string  &dependencies);

  /**
   * \brief Insert dependencies for the value term of the RHS equation of the variable at
   * `index`.
   *
   * \param index Index of variable
   * \param dependencies Container containing list of dependency strings for
   * variable at `index` Hint: {"variable", "grad(variable)", "hess(variable)"}
   */
  template <typename Iterable>
  void
  insert_dependencies_value_term_RHS(const unsigned int &index,
                                     const Iterable     &dependencies);

  /**
   * \brief Insert dependencies for the gradient term of the RHS equation of the variable
   * at `index`.
   *
   * \param index Index of variable
   * \param dependencies Container containing list of dependency strings for
   * variable at `index` Hint: {"variable", "grad(variable)", "hess(variable)"}
   */
  template <typename Iterable>
  void
  insert_dependencies_gradient_term_RHS(const unsigned int &index,
                                        const Iterable     &dependencies);

  /**
   * \brief Insert dependencies for the value term of the LHS equation of the variable at
   * `index`.
   *
   * \param index Index of variable
   * \param dependencies Container containing list of dependency strings for
   * variable at `index` Hint: {"variable", "grad(variable)", "hess(variable)"}
   */
  template <typename Iterable>
  void
  insert_dependencies_value_term_LHS(const unsigned int &index,
                                     const Iterable     &dependencies);

  /**
   * \brief Insert dependencies for the gradient term of the LHS equation of the variable
   * at `index`.
   *
   * \param index Index of variable
   * \param dependencies Container containing list of dependency strings for
   * variable at `index` Hint: {"variable", "grad(variable)", "hess(variable)"}
   */
  template <typename Iterable>
  void
  insert_dependencies_gradient_term_LHS(const unsigned int &index,
                                        const Iterable     &dependencies);

private:
  /**
   * \brief The solutions variable & postprocessing variable attributes
   */
  std::map<unsigned int, variableAttributes> var_attributes;

  /**
   * \brief Perform a suite of assertions on the attributes to ensure that
   * the user's inputs are well-formed
   */
  void
  validate_attributes();

  /**
   * \brief Validate that the variable name is not empty and does not contain any
   * forbidden substrings (names).
   */
  void
  validate_variable_name(const std::string           &name,
                         const std::set<std::string> &forbidden_names,
                         const std::string           &context,
                         unsigned int                 index);

  /**
   * \brief Populate dependencies that we should expect from the user.
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
   * \brief Validate the dependencies (RHS or LHS) that the user has provided.
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
   * \brief Validate the old solution dependencies that the user has provided. There are
   * two criterion: 1. The storage of old fields must be sequential (old_1 must be
   * included before old_2). 2. Old fields cannot be stored for constant equations.
   */
  void
  validate_old_solution_dependencies();

  /**
   * \brief Utility to remove whitespace from strings
   */
  static std::string
  strip_whitespace(const std::string &text);
  // The above function should be moved to a 'utilities' module
};

// Template derived class for variableAttributeLoader for applications.
// `loadVariableAttributes()` and `loadPostProcessorVariableAttributes()` are should be
// filled out in all the applications.
class customAttributeLoader : public variableAttributeLoader
{
public:
  ~customAttributeLoader() override = default;

  void
  loadVariableAttributes() override;
};

PRISMS_PF_END_NAMESPACE

#endif