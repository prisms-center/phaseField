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
 * @brief Class the stores and manages user-defined constants.
 */
template <unsigned int dim>
class UserConstants
{
public:
  using InputVariant = boost::variant<double,
                                      int,
                                      bool,
                                      dealii::Tensor<1, dim>,
                                      dealii::Tensor<2, dim>,
                                      dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>>;

  /**
   * @brief Assign the specified user constant to whatever type.
   */
  InputVariant
  construct_user_constant(std::vector<std::string> &model_constants_strings);

  /**
   * @brief Retrieve the double from the `model_constants` that are defined from the
   * parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * @param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] double
  get_model_constant_double(const std::string &constant_name) const;

  /**
   * @brief Retrieve the int from the `model_constants` that are defined from the
   * parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * @param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] int
  get_model_constant_int(const std::string &constant_name) const;

  /**
   * @brief Retrieve the bool from the `model_constants` that are defined from the
   * parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * @param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] bool
  get_model_constant_bool(const std::string &constant_name) const;

  /**
   * @brief Retrieve the rank 1 tensor from the `model_constants` that are defined from
   * the parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * @param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] dealii::Tensor<1, dim>
  get_model_constant_rank_1_tensor(const std::string &constant_name) const;

  /**
   * @brief Retrieve the rank 2 tensor from the `model_constants` that are defined from
   * the parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * @param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] dealii::Tensor<2, dim>
  get_model_constant_rank_2_tensor(const std::string &constant_name) const;

  /**
   * @brief Retrieve the elasticity tensor from the `model_constants` that are defined
   * from the parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * @param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
  get_model_constant_elasticity_tensor(const std::string &constant_name) const;

  /**
   * @brief Print all user-specified constants
   */
  void
  add_user_constant(const std::string        &constant_name,
                    std::vector<std::string> &model_constants_strings)
  {
    model_constants[constant_name] = construct_user_constant(model_constants_strings);
  };

  /**
   * @brief Print all user-specified constants
   */
  void
  print() const;

private:
  /**
   * @brief Compute the number of tensor rows.
   */
  unsigned int
  compute_tensor_parentheses(const unsigned int             &n_elements,
                             const std::vector<std::string> &tensor_elements);

  /**
   * @brief Remove and leading and trailing parentheses.
   */
  void
  remove_parentheses(std::vector<std::string> &tensor_elements);

  /**
   * @brief Compute a 1st rank tensor from user inputs .
   */
  dealii::Tensor<1, dim>
  compute_rank_1_tensor_constant(const unsigned int             &n_elements,
                                 const std::vector<std::string> &tensor_elements);

  /**
   * @brief Compute a 2nd rank tensor from user inputs .
   */
  dealii::Tensor<2, dim>
  compute_rank_2_tensor_constant(const unsigned int             &n_elements,
                                 const std::vector<std::string> &tensor_elements);

  /**
   * @brief Assign the primitive user constants (e.g., int, double, bool).
   */
  InputVariant
  primitive_model_constant(std::vector<std::string> &model_constants_strings);

  [[nodiscard]] dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
  get_cij_tensor(std::vector<double> elastic_constants,
                 const std::string  &elastic_const_symmetry) const;

  [[nodiscard]] dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
  get_cij_matrix(const ElasticityModel     &model,
                 const std::vector<double> &constants) const;

  /**
   * @brief List of user-defined constants.
   */
  std::map<std::string, InputVariant> model_constants;

  /**
   * @brief Class for printing of variant types. This is bad practice and should be fixed.
   */
  class VariantPrinter : public boost::static_visitor<>
  {
  public:
    void
    operator()(double value) const
    {
      ConditionalOStreams::pout_summary() << value;
    }

    void
    operator()(int value) const
    {
      ConditionalOStreams::pout_summary() << value;
    }

    void
    operator()(bool value) const
    {
      ConditionalOStreams::pout_summary() << std::boolalpha << value;
    }

    void
    operator()(const dealii::Tensor<1, dim> &value) const
    {
      ConditionalOStreams::pout_summary() << "Tensor<1, " << dim << ">: ";
      for (unsigned int i = 0; i < dim; ++i)
        {
          ConditionalOStreams::pout_summary() << value[i] << ' ';
        }
    }

    void
    operator()(const dealii::Tensor<2, dim> &value) const
    {
      ConditionalOStreams::pout_summary() << "Tensor<2, " << dim << ">: ";
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int j = 0; j < dim; ++j)
            {
              ConditionalOStreams::pout_summary() << value[i][j] << ' ';
            }
        }
    }

    template <unsigned int D = dim>
    void
    operator()(const dealii::Tensor<2, (2 * D) - 1 + (D / 3)> &value) const
    requires((D != ((2 * D) - 1 + (D / 3))))
    {
      constexpr unsigned int dimension = (2 * D) - 1 + (D / 3);
      ConditionalOStreams::pout_summary() << "Tensor<2, " << dimension << ">: ";
      for (unsigned int i = 0; i < dimension; ++i)
        {
          for (unsigned int j = 0; j < dimension; ++j)
            {
              ConditionalOStreams::pout_summary() << value[i][j] << ' ';
            }
        }
    }
  };
};

template <unsigned int dim>
inline double
UserConstants<dim>::get_model_constant_double(const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "Mismatch between constants in parameters.prm and CustomPDE.h. The constant "
           "that you attempted to access was " +
           constant_name + "."));

  return boost::get<double>(model_constants.at(constant_name));
}

template <unsigned int dim>
inline int
UserConstants<dim>::get_model_constant_int(const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "Mismatch between constants in parameters.prm and CustomPDE.h. The constant "
           "that you attempted to access was " +
           constant_name + "."));

  return boost::get<int>(model_constants.at(constant_name));
}

template <unsigned int dim>
inline bool
UserConstants<dim>::get_model_constant_bool(const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "Mismatch between constants in parameters.prm and CustomPDE.h. The constant "
           "that you attempted to access was " +
           constant_name + "."));

  return boost::get<bool>(model_constants.at(constant_name));
}

template <unsigned int dim>
inline dealii::Tensor<1, dim>
UserConstants<dim>::get_model_constant_rank_1_tensor(
  const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           " Mismatch between constants in parameters.prm and "
           "CustomPDE.h. The constant that you attempted to access was " +
           constant_name + "."));

  return boost::get<dealii::Tensor<1, dim>>(model_constants.at(constant_name));
}

template <unsigned int dim>
inline dealii::Tensor<2, dim>
UserConstants<dim>::get_model_constant_rank_2_tensor(
  const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "Mismatch between constants in parameters.prm and CustomPDE.h. The constant "
           "that you attempted to access was " +
           constant_name + "."));

  return boost::get<dealii::Tensor<2, dim>>(model_constants.at(constant_name));
}

template <unsigned int dim>
inline dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
UserConstants<dim>::get_model_constant_elasticity_tensor(
  const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "Mismatch between constants in parameters.prm and CustomPDE.h. The constant "
           "that you attempted to access was " +
           constant_name + "."));

  return boost::get<dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>>(
    model_constants.at(constant_name));
}

template <unsigned int dim>
inline unsigned int
UserConstants<dim>::compute_tensor_parentheses(
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

template <unsigned int dim>
inline void
UserConstants<dim>::remove_parentheses(std::vector<std::string> &tensor_elements)
{
  for (std::string &element : tensor_elements)
    {
      boost::range::remove_erase(element, '(');
      boost::range::remove_erase(element, ')');
    }
}

template <unsigned int dim>
inline dealii::Tensor<1, dim>
UserConstants<dim>::compute_rank_1_tensor_constant(
  const unsigned int             &n_elements,
  const std::vector<std::string> &tensor_elements)
{
  AssertThrow(n_elements == 3,
              dealii::ExcMessage("The columns in user-defined constant tensors must be "
                                 "equal to the maximum number of dimensions."));

  dealii::Tensor<1, dim> temp;
  for (unsigned int i = 0; i < dim; i++)
    {
      temp[i] = dealii::Utilities::string_to_double(tensor_elements.at(i));
    }

  return temp;
}

template <unsigned int dim>
inline dealii::Tensor<2, dim>
UserConstants<dim>::compute_rank_2_tensor_constant(
  const unsigned int             &n_elements,
  const std::vector<std::string> &tensor_elements)
{
  AssertThrow(n_elements == 9,
              dealii::ExcMessage("User-defined constant tensor does not have the "
                                 "correct number of elements, matrices must be 3x3."));

  const unsigned int row_length = 3;

  dealii::Tensor<2, dim> temp;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          temp[i][j] =
            dealii::Utilities::string_to_double(tensor_elements.at((i * row_length) + j));
        }
    }

  return temp;
}

template <unsigned int dim>
inline typename UserConstants<dim>::InputVariant
UserConstants<dim>::construct_user_constant(
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
        get_cij_tensor(temp_elastic_constants, elastic_const_symmetry);
      return temp;
    }

  AssertThrow(false,
              dealii::ExcMessage(
                "Only user-defined constant tensors may have multiple elements."));
  return 0;
}

template <unsigned int dim>
inline typename UserConstants<dim>::InputVariant
UserConstants<dim>::primitive_model_constant(
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

template <unsigned int dim>
inline dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
UserConstants<dim>::get_cij_tensor(std::vector<double> elastic_constants,
                                   const std::string  &elastic_const_symmetry) const
{
  // First set the material model
  ElasticityModel mat_model = Isotropic;
  if (elastic_const_symmetry == "isotropic")
    {
      mat_model = ElasticityModel::Isotropic;
    }
  else if (elastic_const_symmetry == "transverse")
    {
      mat_model = ElasticityModel::Transverse;
    }
  else if (elastic_const_symmetry == "orthotropic")
    {
      mat_model = ElasticityModel::Orthotropic;
    }
  else if (elastic_const_symmetry == "anisotropic")
    {
      mat_model = ElasticityModel::Anisotropic;
    }
  else
    {
      AssertThrow(false, dealii::ExcMessage("Invalid elasticity tensor type"));
    }

  // If the material model is anisotropic for a 2D calculation but the elastic
  // constants are given for a 3D calculation, change the elastic constant
  // vector to the 2D form
  constexpr unsigned int max_number = 21;
  if ((mat_model == Anisotropic) && (dim == 2) && elastic_constants.size() == max_number)
    {
      std::vector<double> elastic_constants_temp = elastic_constants;
      elastic_constants.clear();
      const std::vector<unsigned int> indices_2d = {0, 1, 5, 6, 10, 14};
      std::transform(indices_2d.begin(),
                     indices_2d.end(),
                     std::back_inserter(elastic_constants),
                     [&elastic_constants_temp](unsigned int index)
                     {
                       return elastic_constants_temp.at(index);
                     });
    }

  return get_cij_matrix(mat_model, elastic_constants);
}

template <unsigned int dim>
inline dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
UserConstants<dim>::get_cij_matrix(const ElasticityModel     &model,
                                   const std::vector<double> &constants) const
{
  // Initialize compliance tensor
  dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)> compliance;

  switch (dim)
    {
      case 1:
        {
          const int xx_dir = 0;

          switch (model)
            {
                // For the 1D case, it make little sense to accept anything besides an
                // isotropic elasticity tensor. One hiccup, is that if a user is debugging
                // an application and switches to 1D, they will have to modify the
                // elasticity constants to align with that. While more burdensome, there's
                // less of a chance of producing spurious behavior.
              case Isotropic:
                {
                  const double modulus = constants.at(0);

                  compliance[xx_dir][xx_dir] = modulus;
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
          // The voigt indexing scheme for 2 dimensions
          const int xx_dir = 0;
          const int yy_dir = 1;
          const int xy_dir = 2;

          switch (model)
            {
              // Like the 1D case, it is nonsensical to have transverse or orthotropic
              // compliance tensors, so we throw an error.
              case Isotropic:
                {
                  // For isotropic compliance tensors, we can simplify the computation to
                  // two parameters: $\lambda$ and $\mu$, where $\mu$ is the shear
                  // modulus. In cartesian coordinates,
                  // $$$
                  // c_{ijkl} = \lambda \delta_{ij} \delta_{kl} + \mu (\delta_{ik}
                  // \delta_{jl} + \delta_{il} \delta_{kj})
                  // $$$
                  const double modulus = constants.at(0);
                  const double poisson = constants.at(1);

                  const double shear_modulus = modulus / (2 * (1 + poisson));
                  const double lambda =
                    poisson * modulus / ((1 + poisson) * (1 - 2 * poisson));

                  compliance[xx_dir][xx_dir] = compliance[yy_dir][yy_dir] =
                    lambda + 2 * shear_modulus;
                  compliance[xy_dir][xy_dir] = shear_modulus;
                  compliance[xx_dir][yy_dir] = compliance[yy_dir][xx_dir] = lambda;
                  break;
                }
              case Anisotropic:
                {
                  // In the anisotropic case, every entry is specified (given the symmetry
                  // constraints). Also, ignore magic numbers because it is simpler to
                  // hardcode this.

                  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)

                  compliance[xx_dir][xx_dir] = constants.at(0);
                  compliance[yy_dir][yy_dir] = constants.at(1);
                  compliance[xy_dir][xy_dir] = constants.at(2);
                  compliance[xx_dir][yy_dir] = compliance[yy_dir][xx_dir] =
                    constants.at(3);
                  compliance[xx_dir][xy_dir] = compliance[xy_dir][xx_dir] =
                    constants.at(4);
                  compliance[yy_dir][xy_dir] = compliance[xy_dir][yy_dir] =
                    constants.at(5);

                  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
                  break;
                }
              default:
                AssertThrow(false, dealii::ExcMessage("Invalid elasticity model type"));
            }
          break;
        }
      case 3:
        {
          const int xx_dir = 0;
          const int yy_dir = 1;
          const int zz_dir = 2;
          const int yz_dir = 3;
          const int xz_dir = 4;
          const int xy_dir = 5;

          switch (model)
            {
              case Isotropic:
                {
                  // For isotropic compliance tensors, we can simplify the computation to
                  // two parameters: $\lambda$ and $\mu$, where $\mu$ is the shear
                  // modulus. In cartesian coordinates,
                  // $$$
                  // c_{ijkl} = \lambda \delta_{ij} \delta_{kl} + \mu (\delta_{ik}
                  // \delta_{jl} + \delta_{il} \delta_{kj})
                  // $$$
                  const double modulus = constants.at(0);
                  const double poisson = constants.at(1);

                  const double shear_modulus = modulus / (2 * (1 + poisson));
                  const double lambda =
                    poisson * modulus / ((1 + poisson) * (1 - 2 * poisson));

                  compliance[xx_dir][xx_dir]     = compliance[yy_dir][yy_dir] =
                    compliance[zz_dir][zz_dir]   = lambda + 2 * shear_modulus;
                  compliance[yz_dir][yz_dir]     = compliance[xz_dir][xz_dir] =
                    compliance[xy_dir][xy_dir]   = shear_modulus;
                  compliance[xx_dir][yy_dir]     = compliance[yy_dir][xx_dir] =
                    compliance[xx_dir][zz_dir]   = compliance[zz_dir][xx_dir] =
                      compliance[yy_dir][zz_dir] = compliance[zz_dir][yy_dir] = lambda;
                  break;
                }
              case Anisotropic:
                {
                  // In the anisotropic case, every entry is specified (given the symmetry
                  // constraints). Also, ignore magic numbers because it is simpler to
                  // hardcode this.

                  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)

                  compliance[xx_dir][xx_dir] = constants[0];
                  compliance[yy_dir][yy_dir] = constants[1];
                  compliance[zz_dir][zz_dir] = constants[2];
                  compliance[yz_dir][yz_dir] = constants[3];
                  compliance[xz_dir][xz_dir] = constants[4];
                  compliance[xy_dir][xy_dir] = constants[5];
                  compliance[xx_dir][yy_dir] = compliance[yy_dir][xx_dir] = constants[6];
                  compliance[xx_dir][zz_dir] = compliance[zz_dir][xx_dir] = constants[7];
                  compliance[xx_dir][yz_dir] = compliance[yz_dir][xx_dir] = constants[8];
                  compliance[xx_dir][xz_dir] = compliance[xz_dir][xx_dir] = constants[9];
                  compliance[xx_dir][xy_dir] = compliance[xy_dir][xx_dir] = constants[10];
                  compliance[yy_dir][zz_dir] = compliance[zz_dir][yy_dir] = constants[11];
                  compliance[yy_dir][yz_dir] = compliance[yz_dir][yy_dir] = constants[12];
                  compliance[yy_dir][xz_dir] = compliance[xz_dir][yy_dir] = constants[13];
                  compliance[yy_dir][xy_dir] = compliance[xy_dir][yy_dir] = constants[14];
                  compliance[zz_dir][yz_dir] = compliance[yz_dir][zz_dir] = constants[15];
                  compliance[zz_dir][xz_dir] = compliance[xz_dir][zz_dir] = constants[16];
                  compliance[zz_dir][xy_dir] = compliance[xy_dir][zz_dir] = constants[17];
                  compliance[yz_dir][xz_dir] = compliance[xz_dir][yz_dir] = constants[18];
                  compliance[yz_dir][xy_dir] = compliance[xy_dir][yz_dir] = constants[19];
                  compliance[xz_dir][xy_dir] = compliance[xy_dir][xz_dir] = constants[20];

                  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
                  break;
                }
              case Transverse:
              // TODO (landinjm): implement
              case Orthotropic:
              // TODO (landinjm): implement
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

  return compliance;
}

template <unsigned int dim>
void
UserConstants<dim>::print() const
{
  if (!model_constants.empty())
    {
      ConditionalOStreams::pout_summary()
        << "================================================\n"
        << "  User Constants\n"
        << "================================================\n";

      for (const auto &[constant_name, variant] : model_constants)
        {
          ConditionalOStreams::pout_summary() << constant_name << ": ";
          boost::apply_visitor(VariantPrinter(), variant);
          ConditionalOStreams::pout_summary() << "\n";
        }
      ConditionalOStreams::pout_summary() << "\n" << std::flush;
    }
}

PRISMS_PF_END_NAMESPACE
