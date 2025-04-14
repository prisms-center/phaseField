// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/variant.hpp>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

#include <map>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Class the stores and manages user-defined constants.
 */
template <int dim>
class userConstants
{
public:
  using InputVariant = boost::variant<double,
                                      int,
                                      bool,
                                      dealii::Tensor<1, dim>,
                                      dealii::Tensor<2, dim>,
                                      dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>>;

  /**
   * \brief Assign the specified user constant to whatever type.
   */
  InputVariant
  construct_user_constant(std::vector<std::string> &model_constants_strings);

  /**
   * \brief Retrieve the double from the `model_constants` that are defined from the
   * parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] double
  get_model_constant_double(const std::string &constant_name) const;

  /**
   * \brief Retrieve the int from the `model_constants` that are defined from the
   * parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] int
  get_model_constant_int(const std::string &constant_name) const;

  /**
   * \brief Retrieve the bool from the `model_constants` that are defined from the
   * parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] bool
  get_model_constant_bool(const std::string &constant_name) const;

  /**
   * \brief Retrieve the rank 1 tensor from the `model_constants` that are defined from
   * the parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] dealii::Tensor<1, dim>
  get_model_constant_rank_1_tensor(const std::string &constant_name) const;

  /**
   * \brief Retrieve the rank 2 tensor from the `model_constants` that are defined from
   * the parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] dealii::Tensor<2, dim>
  get_model_constant_rank_2_tensor(const std::string &constant_name) const;

  /**
   * \brief Retrieve the elasticity tensor from the `model_constants` that are defined
   * from the parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
  get_model_constant_elasticity_tensor(const std::string &constant_name) const;

  // List of user-defined constants
  std::map<std::string, InputVariant> model_constants;

private:
  /**
   * \brief Compute the number of tensor rows.
   */
  unsigned int
  compute_tensor_parentheses(const unsigned int             &n_elements,
                             const std::vector<std::string> &tensor_elements);

  /**
   * \brief Remove and leading and trailing parentheses.
   */
  void
  remove_parentheses(std::vector<std::string> &tensor_elements);

  /**
   * \brief Compute a 1st rank tensor from user inputs .
   */
  dealii::Tensor<1, dim>
  compute_rank_1_tensor_constant(const unsigned int             &n_elements,
                                 const std::vector<std::string> &tensor_elements);

  /**
   * \brief Compute a 2nd rank tensor from user inputs .
   */
  dealii::Tensor<2, dim>
  compute_rank_2_tensor_constant(const unsigned int             &n_elements,
                                 const std::vector<std::string> &tensor_elements);

  /**
   * \brief Assign the primitive user constants (e.g., int, double, bool).
   */
  InputVariant
  primitive_model_constant(std::vector<std::string> &model_constants_strings);

  [[nodiscard]] dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
  get_Cij_tensor(std::vector<double> elastic_constants,
                 const std::string  &elastic_const_symmetry) const;

  [[nodiscard]] dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
  getCIJMatrix(const elasticityModel &model, const std::vector<double> &constants) const;
};

template <int dim>
inline double
userConstants<dim>::get_model_constant_double(const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "Mismatch between constants in parameters.prm and customPDE.h. The constant "
           "that you attempted to access was " +
           constant_name + "."));

  return boost::get<double>(model_constants.at(constant_name));
}

template <int dim>
inline int
userConstants<dim>::get_model_constant_int(const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "Mismatch between constants in parameters.prm and customPDE.h. The constant "
           "that you attempted to access was " +
           constant_name + "."));

  return boost::get<int>(model_constants.at(constant_name));
}

template <int dim>
inline bool
userConstants<dim>::get_model_constant_bool(const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "Mismatch between constants in parameters.prm and customPDE.h. The constant "
           "that you attempted to access was " +
           constant_name + "."));

  return boost::get<bool>(model_constants.at(constant_name));
}

template <int dim>
inline dealii::Tensor<1, dim>
userConstants<dim>::get_model_constant_rank_1_tensor(
  const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           " Mismatch between constants in parameters.prm and "
           "customPDE.h. The constant that you attempted to access was " +
           constant_name + "."));

  return boost::get<dealii::Tensor<1, dim>>(model_constants.at(constant_name));
}

template <int dim>
inline dealii::Tensor<2, dim>
userConstants<dim>::get_model_constant_rank_2_tensor(
  const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "Mismatch between constants in parameters.prm and customPDE.h. The constant "
           "that you attempted to access was " +
           constant_name + "."));

  return boost::get<dealii::Tensor<2, dim>>(model_constants.at(constant_name));
}

template <int dim>
inline dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
userConstants<dim>::get_model_constant_elasticity_tensor(
  const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "Mismatch between constants in parameters.prm and customPDE.h. The constant "
           "that you attempted to access was " +
           constant_name + "."));

  return boost::get<dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>>(
    model_constants.at(constant_name));
}

template <int dim>
inline unsigned int
userConstants<dim>::compute_tensor_parentheses(
  const unsigned int             &n_elements,
  const std::vector<std::string> &tensor_elements)
{
  unsigned int open_parentheses  = 0;
  unsigned int close_parentheses = 0;

  for (unsigned int element = 0; element < n_elements; element++)
    {
      for (const char character : tensor_elements.at(element))
        {
          if (character == '(')
            {
              ++open_parentheses;
            }
          else if (character == ')')
            {
              ++close_parentheses;
            }
        }
    }

  if (open_parentheses != close_parentheses)
    {
      AssertThrow(false,
                  dealii::ExcMessage("User-defined elastic constant list does not have "
                                     "the same number of open and close parentheses."));
    }

  return open_parentheses;
}

template <int dim>
inline void
userConstants<dim>::remove_parentheses(std::vector<std::string> &tensor_elements)
{
  for (std::string &element : tensor_elements)
    {
      boost::range::remove_erase(element, '(');
      boost::range::remove_erase(element, ')');
    }
}

template <int dim>
inline dealii::Tensor<1, dim>
userConstants<dim>::compute_rank_1_tensor_constant(
  const unsigned int             &n_elements,
  const std::vector<std::string> &tensor_elements)
{
  AssertThrow(n_elements == 3,
              dealii::ExcMessage("The columns in user-defined constant tensors must be "
                                 "equal to the maximum number of dimensions."));

  dealii::Tensor<1, dim> temp;
  for (int i = 0; i < dim; i++)
    {
      temp[i] = dealii::Utilities::string_to_double(tensor_elements.at(i));
    }

  return temp;
}

template <int dim>
inline dealii::Tensor<2, dim>
userConstants<dim>::compute_rank_2_tensor_constant(
  const unsigned int             &n_elements,
  const std::vector<std::string> &tensor_elements)
{
  AssertThrow(n_elements == 9,
              dealii::ExcMessage("User-defined constant tensor does not have the "
                                 "correct number of elements, matrices must be 3x3."));

  const unsigned int row_length = dim;

  dealii::Tensor<2, dim> temp;
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          temp[i][j] =
            dealii::Utilities::string_to_double(tensor_elements.at((i * row_length) + j));
        }
    }

  return temp;
}

template <int dim>
inline typename userConstants<dim>::InputVariant
userConstants<dim>::construct_user_constant(
  std::vector<std::string> &model_constants_strings)
{
  // Ensure that the input includes a value and a type
  AssertThrow(
    model_constants_strings.size() > 1,
    dealii::ExcMessage(
      "At least two fields are required for user-defined variables (value and type)."));

  std::vector<std::string> model_constants_type_strings =
    dealii::Utilities::split_string_list(model_constants_strings.at(
                                           model_constants_strings.size() - 1),
                                         ' ');

  if (model_constants_strings.size() == 2)
    {
      return primitive_model_constant(model_constants_strings);
    }

  if (boost::iequals(model_constants_type_strings.at(0), "tensor"))
    {
      const unsigned int n_elements = model_constants_strings.size() - 1;

      const unsigned int open_parentheses =
        compute_tensor_parentheses(n_elements, model_constants_strings);

      AssertThrow(open_parentheses <= 4,
                  FeatureNotImplemented("3rd rank tensors and above"));

      remove_parentheses(model_constants_strings);

      // Rank 1 tensor
      if (open_parentheses == 1)
        {
          return compute_rank_1_tensor_constant(n_elements, model_constants_strings);
        }
      // Rank 2 tensor
      return compute_rank_2_tensor_constant(n_elements, model_constants_strings);
    }
  if (boost::iequals(model_constants_type_strings.at(1), "elastic") &&
      boost::iequals(model_constants_type_strings.at(2), "constants"))
    {
      const unsigned int n_elements = model_constants_strings.size() - 1;

      remove_parentheses(model_constants_strings);

      // Load in the elastic constants as a vector
      std::vector<double> temp_elastic_constants(n_elements);
      for (unsigned int i = 0; i < n_elements; i++)
        {
          temp_elastic_constants[i] =
            dealii::Utilities::string_to_double(model_constants_strings.at(i));
        }
      const std::string &elastic_const_symmetry = model_constants_type_strings.at(0);
      dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)> temp =
        get_Cij_tensor(temp_elastic_constants, elastic_const_symmetry);
      return temp;
    }

  AssertThrow(false,
              dealii::ExcMessage(
                "Only user-defined constant tensors may have multiple elements."));
  return 0;
}

template <int dim>
inline typename userConstants<dim>::InputVariant
userConstants<dim>::primitive_model_constant(
  std::vector<std::string> &model_constants_strings)
{
  std::vector<std::string> model_constants_type_strings =
    dealii::Utilities::split_string_list(model_constants_strings.at(
                                           model_constants_strings.size() - 1),
                                         ' ');

  if (boost::iequals(model_constants_type_strings.at(0), "double"))
    {
      return dealii::Utilities::string_to_double(model_constants_strings.at(0));
    }
  if (boost::iequals(model_constants_type_strings.at(0), "int"))
    {
      return dealii::Utilities::string_to_int(model_constants_strings.at(0));
    }
  if (boost::iequals(model_constants_type_strings.at(0), "bool"))
    {
      return boost::iequals(model_constants_strings.at(0), "true");
    }

  AssertThrow(false,
              dealii::ExcMessage(
                "The type for user-defined variables must be `double`, `int`, "
                "`bool`, `tensor`, or `elastic constants`."));
  return 0;
}

template <int dim>
inline dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
userConstants<dim>::get_Cij_tensor(std::vector<double> elastic_constants,
                                   const std::string  &elastic_const_symmetry) const
{
  // First set the material model
  elasticityModel mat_model = ISOTROPIC;
  if (elastic_const_symmetry == "isotropic")
    {
      mat_model = elasticityModel::ISOTROPIC;
    }
  else if (elastic_const_symmetry == "transverse")
    {
      mat_model = elasticityModel::TRANSVERSE;
    }
  else if (elastic_const_symmetry == "orthotropic")
    {
      mat_model = elasticityModel::ORTHOTROPIC;
    }
  else if (elastic_const_symmetry == "anisotropic")
    {
      mat_model = elasticityModel::ANISOTROPIC;
    }
  else
    {
      AssertThrow(false, dealii::ExcMessage("Invalid elasticity tensor type"));
    }

  // If the material model is anisotropic for a 2D calculation but the elastic
  // constants are given for a 3D calculation, change the elastic constant
  // vector to the 2D form
  constexpr unsigned int max_number = 21;
  if ((mat_model == ANISOTROPIC) && (dim == 2) && elastic_constants.size() == max_number)
    {
      std::vector<double> elastic_constants_temp = elastic_constants;
      elastic_constants.clear();
      const std::vector<unsigned int> indices_2D = {0, 1, 5, 6, 10, 14};
      std::transform(indices_2D.begin(),
                     indices_2D.end(),
                     std::back_inserter(elastic_constants),
                     [&elastic_constants_temp](unsigned int index)
                     {
                       return elastic_constants_temp.at(index);
                     });
    }

  return getCIJMatrix(mat_model, elastic_constants);
}

template <int dim>
inline dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
userConstants<dim>::getCIJMatrix(const elasticityModel     &model,
                                 const std::vector<double> &constants) const
{
  // Initialize tensor
  dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)> CIJ;

  switch (dim)
    {
      case 1:
        {
          switch (model)
            {
              // TODO (landinjm): Should we both fixing this for the other cases and just
              // selecting the x index. It would allow the user to switch from 2D to 1D to
              // 3D, but would produce unexpected behavior.
              case ISOTROPIC:
                {
                  const double modulus = constants.at(0);

                  CIJ[0][0] = modulus;
                  break;
                }
              default:
                AssertThrow(false,
                            dealii::ExcMessage(
                              "Invalid elasticity model type for 1D. We only accept "
                              "isotropic elasticity tensors."));
            }
          break;
        }
      case 2:
        {
          switch (model)
            {
              case ISOTROPIC:
                {
                  // TODO (landinjm): Document this
                  const double modulus = constants.at(0);
                  const double poisson = constants.at(1);

                  const double shear_modulus = modulus / (2 * (1 + poisson));
                  const double lambda =
                    poisson * modulus / ((1 + poisson) * (1 - 2 * poisson));

                  CIJ[0][0] = lambda + 2 * shear_modulus;
                  CIJ[1][1] = lambda + 2 * shear_modulus;
                  CIJ[2][2] = shear_modulus;
                  CIJ[0][1] = lambda;
                  CIJ[1][0] = lambda;
                  break;
                }
              case ANISOTROPIC:
              default:
                AssertThrow(false, dealii::ExcMessage("Invalid elasticity model type"));
            }
          break;
        }
      case 3:
        {
          switch (model)
            {
              case ISOTROPIC:
              case TRANSVERSE:
              case ORTHOTROPIC:
              case ANISOTROPIC:
              default:
                AssertThrow(false, dealii::ExcMessage("Invalid elasticity model type"));
            }
          break;
        }
      default:
        {
          Assert(false, UnreachableCode());
        }
    }

  return CIJ;
}

PRISMS_PF_END_NAMESPACE
