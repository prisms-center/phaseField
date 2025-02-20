// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef user_constants_h
#define user_constants_h

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/variant.hpp>

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>

#include <iomanip>
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

  dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
  getCIJMatrix(const elasticityModel &model, const std::vector<double> &constants) const;
};

template <int dim>
inline double
userConstants<dim>::get_model_constant_double(const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
           "customPDE.h. The constant that you attempted to access was " +
           constant_name + "."));

  return boost::get<double>(model_constants.at(constant_name));
}

template <int dim>
inline int
userConstants<dim>::get_model_constant_int(const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
           "customPDE.h. The constant that you attempted to access was " +
           constant_name + "."));

  return boost::get<int>(model_constants.at(constant_name));
}

template <int dim>
inline bool
userConstants<dim>::get_model_constant_bool(const std::string &constant_name) const
{
  Assert(model_constants.find(constant_name) != model_constants.end(),
         dealii::ExcMessage(
           "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
           "customPDE.h. The constant that you attempted to access was " +
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
           "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
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
           "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
           "customPDE.h. The constant that you attempted to access was " +
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
           "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
           "customPDE.h. The constant that you attempted to access was " +
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
      for (const char c : tensor_elements.at(element))
        {
          if (c == '(')
            {
              ++open_parentheses;
            }
          else if (c == ')')
            {
              ++close_parentheses;
            }
        }
    }

  if (open_parentheses != close_parentheses)
    {
      std::cerr << "PRISMS-PF ERROR: User-defined elastic constant "
                   "list does not have the same number of open and "
                   "close parentheses.\n";
      abort();
    }

  return open_parentheses;
}

template <int dim>
inline void
userConstants<dim>::remove_parentheses(std::vector<std::string> &tensor_elements)
{
  for (std::string &element : tensor_elements)
    {
      element.erase(std::remove(element.begin(), element.end(), '('), element.end());
      element.erase(std::remove(element.begin(), element.end(), ')'), element.end());
    }
}

template <int dim>
inline dealii::Tensor<1, dim>
userConstants<dim>::compute_rank_1_tensor_constant(
  const unsigned int             &n_elements,
  const std::vector<std::string> &tensor_elements)
{
  AssertThrow(n_elements > 1 && n_elements < 4,
              dealii::ExcMessage("PRISMS-PF Error: The columns in user-defined constant "
                                 "tensors cannot be longer than 3 elements (internally "
                                 "truncated to the number of dimensions)."));

  dealii::Tensor<1, dim> temp;
  for (unsigned int i = 0; i < dim; i++)
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
  unsigned int row_length = 0;
  if (n_elements == 4)
    {
      AssertThrow(dim < 3,
                  dealii::ExcMessage(
                    "PRISMS-PF ERROR: User-defined constant tensor does not have "
                    "enough elements. For 3D calculations matrices must be 3x3."));

      row_length = 2;
    }
  else if (n_elements == 9)
    {
      row_length = 3;
    }
  else
    {
      AssertThrow(false,
                  dealii::ExcMessage("PRISMS-PF ERROR: User-defined constant tensor does "
                                     "not have the correct number of elements, matrices "
                                     "must be 2x2 or 3x3."));
    }

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

template <int dim>
inline typename userConstants<dim>::InputVariant
userConstants<dim>::construct_user_constant(
  std::vector<std::string> &model_constants_strings)
{
  // Ensure that the input includes a value and a type
  AssertThrow(model_constants_strings.size() > 1,
              dealii::ExcMessage("PRISMS-PF Error: At least two fields are required for "
                                 "user-defined variables (value and type)."));

  std::vector<std::string> model_constants_type_strings =
    dealii::Utilities::split_string_list(model_constants_strings.at(
                                           model_constants_strings.size() - 1),
                                         ' ');

  if (model_constants_strings.size() == 2)
    {
      return primitive_model_constant(model_constants_strings);
    }
  else
    {
      if (boost::iequals(model_constants_type_strings.at(0), "tensor"))
        {
          const unsigned int n_elements = model_constants_strings.size() - 1;

          const unsigned int open_parentheses =
            compute_tensor_parentheses(n_elements, model_constants_strings);
          remove_parentheses(model_constants_strings);

          // Rank 1 tensor
          if (open_parentheses < 3)
            {
              return compute_rank_1_tensor_constant(n_elements, model_constants_strings);
            }
          // Rank 2 tensor
          else if (open_parentheses < 5)
            {
              return compute_rank_2_tensor_constant(n_elements, model_constants_strings);
            }
        }
      else if (boost::iequals(model_constants_type_strings.at(1), "elastic") &&
               boost::iequals(model_constants_type_strings.at(2), "constants"))
        {
          const unsigned int n_elements = model_constants_strings.size() - 1;

          remove_parentheses(model_constants_strings);

          // Load in the elastic constants as a vector
          std::vector<double> temp_elastic_constants;
          for (unsigned int i = 0; i < n_elements; i++)
            {
              temp_elastic_constants.push_back(
                dealii::Utilities::string_to_double(model_constants_strings.at(i)));
            }

          const std::string &elastic_const_symmetry = model_constants_type_strings.at(0);
          dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)> temp =
            get_Cij_tensor(temp_elastic_constants, elastic_const_symmetry);
          return temp;
        }
      else
        {
          AssertThrow(false,
                      dealii::ExcMessage(
                        "PRISMS-PF ERROR: Only user-defined constant tensors may "
                        "have multiple elements."));
        }
      return 0;
    }
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
  else if (boost::iequals(model_constants_type_strings.at(0), "int"))
    {
      return dealii::Utilities::string_to_int(model_constants_strings.at(0));
    }
  else if (boost::iequals(model_constants_type_strings.at(0), "bool"))
    {
      bool temp = boost::iequals(model_constants_strings.at(0), "true");
      return temp;
    }
  else
    {
      AssertThrow(false,
                  dealii::ExcMessage(
                    "PRISMS-PF Error: The type for user-defined variables must be "
                    "`double`, `int`, `bool`, `tensor`, or `elastic constants`."));
      return 0;
    }
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
      // Should change to an exception
      std::cerr << "Elastic material model is invalid, please use isotropic, "
                   "transverse, orthotropic, or anisotropic\n";
    }

  // If the material model is anisotropic for a 2D calculation but the elastic
  // constants are given for a 3D calculation, change the elastic constant
  // vector to the 2D form
  if ((mat_model == ANISOTROPIC) && (dim == 2) && elastic_constants.size() == 21)
    {
      std::vector<double> elastic_constants_temp = elastic_constants;
      elastic_constants.clear();
      const std::vector<unsigned int> indices_2D = {0, 1, 5, 6, 10, 14};
      for (const auto &index : indices_2D)
        {
          elastic_constants.push_back(elastic_constants_temp.at(index));
        }
    }

  return getCIJMatrix(mat_model, elastic_constants);
}

template <int dim>
inline dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)>
userConstants<dim>::getCIJMatrix(const elasticityModel     &model,
                                 const std::vector<double> &constants) const
{
  // CIJ.fill(0.0);
  dealii::Tensor<2, (2 * dim) - 1 + (dim / 3)> CIJ;

  conditionalOStreams::pout_base() << "Reading material model:";
  switch (dim)
    {
      case 1:
        {
          conditionalOStreams::pout_base() << " 1D ";
          // 1D models
          switch (model)
            {
              case ISOTROPIC:
                {
                  conditionalOStreams::pout_base() << " ISOTROPIC \n";
                  CIJ[0][0] = constants[0];
                  break;
                }
              default:
                {
                  exit(-1);
                }
            }
          break;
        }
      case 2:
        {
          conditionalOStreams::pout_base() << " 2D ";
          // 2D models
          switch (model)
            {
              case ISOTROPIC:
                {
                  conditionalOStreams::pout_base() << " ISOTROPIC \n";
                  const double E      = constants[0];
                  const double nu     = constants[1];
                  const double mu     = E / (2 * (1 + nu));
                  const double lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
                  CIJ[0][0]           = lambda + 2 * mu;
                  CIJ[1][1]           = lambda + 2 * mu;
                  CIJ[2][2]           = mu;
                  CIJ[0][1] = CIJ[1][0] = lambda;
                  break;
                }
              case ANISOTROPIC:
                {
                  conditionalOStreams::pout_base() << " ANISOTROPIC \n";
                  CIJ[0][0] = constants[0];             // C11
                  CIJ[1][1] = constants[1];             // C22
                  CIJ[2][2] = constants[2];             // C33
                  CIJ[0][1] = CIJ[1][0] = constants[3]; // C12
                  CIJ[0][2] = CIJ[2][0] = constants[4]; // C13
                  CIJ[1][2] = CIJ[2][1] = constants[5]; // C23
                  break;
                }
              default:
                {
                  exit(-1);
                }
            }
          break;
        }
      case 3:
        {
          conditionalOStreams::pout_base() << " 3D ";
          // 3D models
          switch (model)
            {
              case ISOTROPIC:
                {
                  conditionalOStreams::pout_base() << " ISOTROPIC \n";
                  const double E      = constants[0];
                  const double nu     = constants[1];
                  const double mu     = E / (2 * (1 + nu));
                  const double lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
                  CIJ[0][0]           = lambda + 2 * mu;
                  CIJ[1][1]           = lambda + 2 * mu;
                  CIJ[2][2]           = lambda + 2 * mu;
                  CIJ[3][3]           = mu;
                  CIJ[4][4]           = mu;
                  CIJ[5][5]           = mu;
                  CIJ[0][1] = CIJ[1][0] = lambda;
                  CIJ[0][2] = CIJ[2][0] = lambda;
                  CIJ[1][2] = CIJ[2][1] = lambda;
                  break;
                }
              case TRANSVERSE:
                {
                  conditionalOStreams::pout_base() << " TRANSVERSE \n";
                  CIJ[0][0] = constants[0];                        // C11
                  CIJ[1][1] = constants[0];                        // C11
                  CIJ[2][2] = constants[1];                        // C33
                  CIJ[3][3] = constants[2];                        // C44
                  CIJ[4][4] = constants[2];                        // C44
                  CIJ[5][5] = (constants[0] - constants[3]) / 2.0; //(C11-C12)/2
                  CIJ[0][1] = CIJ[1][0] = constants[3];            // C12
                  CIJ[0][2] = CIJ[2][0] = constants[4];            // C13
                  CIJ[1][2] = CIJ[2][1] = constants[4];            // C13
                  break;
                }
              case ORTHOTROPIC:
                {
                  conditionalOStreams::pout_base() << " ORTHOTROPIC \n";
                  CIJ[0][0] = constants[0];             // C11
                  CIJ[1][1] = constants[1];             // C22
                  CIJ[2][2] = constants[2];             // C33
                  CIJ[3][3] = constants[3];             // C44
                  CIJ[4][4] = constants[4];             // C55
                  CIJ[5][5] = constants[5];             // C66
                  CIJ[0][1] = CIJ[1][0] = constants[6]; // C12
                  CIJ[0][2] = CIJ[2][0] = constants[7]; // C13
                  CIJ[1][2] = CIJ[2][1] = constants[8]; // C23
                  break;
                }
              case ANISOTROPIC:
                {
                  conditionalOStreams::pout_base() << " ANISOTROPIC \n";
                  CIJ[0][0] = constants[0];              // C11
                  CIJ[1][1] = constants[1];              // C22
                  CIJ[2][2] = constants[2];              // C33
                  CIJ[3][3] = constants[3];              // C44
                  CIJ[4][4] = constants[4];              // C55
                  CIJ[5][5] = constants[5];              // C66
                  CIJ[0][1] = CIJ[1][0] = constants[6];  // C12
                  CIJ[0][2] = CIJ[2][0] = constants[7];  // C13
                  CIJ[0][3] = CIJ[3][0] = constants[8];  // C14
                  CIJ[0][4] = CIJ[4][0] = constants[9];  // C15
                  CIJ[0][5] = CIJ[5][0] = constants[10]; // C16
                  CIJ[1][2] = CIJ[2][1] = constants[11]; // C23
                  CIJ[1][3] = CIJ[3][1] = constants[12]; // C24
                  CIJ[1][4] = CIJ[4][1] = constants[13]; // C25
                  CIJ[1][5] = CIJ[5][1] = constants[14]; // C26
                  CIJ[2][3] = CIJ[3][2] = constants[15]; // C34
                  CIJ[2][4] = CIJ[4][2] = constants[16]; // C35
                  CIJ[2][5] = CIJ[5][2] = constants[17]; // C36
                  CIJ[3][4] = CIJ[4][3] = constants[18]; // C45
                  CIJ[3][5] = CIJ[5][3] = constants[19]; // C46
                  CIJ[4][5] = CIJ[5][4] = constants[20]; // C56
                  break;
                }
              default:
                {
                  exit(-1);
                }
            }
          break;
        }
      default:
        {
          exit(-1);
        }
    }
  // print CIJ to terminal
  conditionalOStreams::pout_base() << "Elasticity matrix (Voigt notation):\n";
  constexpr unsigned int voight_matrix_size = (2 * dim) - 1 + (dim / 3);
  for (unsigned int i = 0; i < voight_matrix_size; i++)
    {
      for (unsigned int j = 0; j < voight_matrix_size; j++)
        {
          conditionalOStreams::pout_base() << std::setw(8) << std::setprecision(3)
                                           << std::scientific << CIJ[i][j] << " ";
        }
      conditionalOStreams::pout_base() << "\n";
    }
  conditionalOStreams::pout_base() << "\n";
  return CIJ;
}

PRISMS_PF_END_NAMESPACE

#endif