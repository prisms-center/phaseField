#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include "../../include/userInputParameters.h"

// ==========================================================================================
// Method to extract the user-defined model constants
// ==========================================================================================
template <int dim>
void
userInputParameters<dim>::load_user_constants(inputFileReader          &input_file_reader,
                                              dealii::ParameterHandler &parameter_handler)
{
  unsigned int number_of_constants = input_file_reader.num_constants;

  for (unsigned int i = 0; i < input_file_reader.model_constant_names.size(); i++)
    {
      model_constant_name_map[input_file_reader.model_constant_names[i]] = i;
    }

  for (unsigned int i = 0; i < number_of_constants; i++)
    {
      std::string constants_text = "Model constant ";
      constants_text.append(input_file_reader.model_constant_names[i]);
      // std::cout << input_file_reader.model_constant_names[i] << std::endl;
      std::vector<std::string> model_constants_strings =
        dealii::Utilities::split_string_list(parameter_handler.get(constants_text));

      // Ensure that the input includes a value and a type
      if (model_constants_strings.size() < 2)
        {
          std::cerr << "PRISMS-PF ERROR: Users must input two fields for "
                       "user-defined variables (value and type)."
                    << std::endl;
          abort();
        }

      std::vector<std::string> model_constants_type_strings =
        dealii::Utilities::split_string_list(model_constants_strings.at(
                                               model_constants_strings.size() - 1),
                                             ' ');

      if (model_constants_strings.size() == 2)
        {
          if (boost::iequals(model_constants_type_strings.at(0), "double"))
            {
              model_constants.push_back(
                dealii::Utilities::string_to_double(model_constants_strings.at(0)));
            }
          else if (boost::iequals(model_constants_type_strings.at(0), "int"))
            {
              model_constants.push_back(
                dealii::Utilities::string_to_int(model_constants_strings.at(0)));
            }
          else if (boost::iequals(model_constants_type_strings.at(0), "bool"))
            {
              bool temp;
              if (boost::iequals(model_constants_strings.at(0), "true"))
                {
                  temp = true;
                }
              else
                {
                  temp = false;
                }
              model_constants.push_back(temp);
            }
          else
            {
              std::cerr << "PRISMS-PF ERROR: The type for user-defined "
                           "variables must be 'double', 'int', 'bool', "
                           "'tensor', or 'elastic constants'."
                        << std::endl;
              abort();
            }
        }
      else
        {
          if (boost::iequals(model_constants_type_strings.at(0), "tensor"))
            {
              unsigned int num_elements = model_constants_strings.size() - 1;

              // Strip parentheses from the input, counting how many rows there
              // are
              unsigned int open_parentheses  = 0;
              unsigned int close_parentheses = 0;
              for (unsigned int element = 0; element < num_elements; element++)
                {
                  std::size_t index = 0;
                  while (index != std::string::npos)
                    {
                      index = model_constants_strings.at(element).find("(");
                      if (index != std::string::npos)
                        {
                          model_constants_strings.at(element).erase(index, 1);
                          open_parentheses++;
                        }
                    }
                  index = 0;
                  while (index != std::string::npos)
                    {
                      index = model_constants_strings.at(element).find(")");
                      if (index != std::string::npos)
                        {
                          model_constants_strings.at(element).erase(index, 1);
                          close_parentheses++;
                        }
                    }
                }
              if (open_parentheses != close_parentheses)
                {
                  std::cerr << "PRISMS-PF ERROR: User-defined constant tensor does not "
                               "have the same number of open and close parentheses."
                            << std::endl;
                  abort();
                }
              // Rank 1 tensor
              if (open_parentheses < 3)
                {
                  if (num_elements > 1 && num_elements < 4)
                    {
                      dealii::Tensor<1, dim> temp;
                      for (unsigned int i = 0; i < dim; i++)
                        {
                          temp[i] = dealii::Utilities::string_to_double(
                            model_constants_strings.at(i));
                        }
                      model_constants.push_back(temp);
                    }
                  else
                    {
                      std::cerr << "PRISMS-PF ERROR: The columns in user-defined "
                                   "constant tensors cannot be longer than 3 elements "
                                   "(internally truncated to the number of dimensions)."
                                << std::endl;
                      abort();
                    }
                }
              // Rank 2 tensor
              else if (open_parentheses < 5)
                {
                  unsigned int row_length;
                  if (num_elements == 4)
                    {
                      row_length = 2;
                      if (dim > 2)
                        {
                          std::cerr << "PRISMS-PF ERROR: User-defined constant tensor "
                                       "does not have enough elements, for 3D "
                                       "calculations matrices must be 3x3."
                                    << std::endl;
                          abort();
                        }
                    }
                  else if (num_elements == 9)
                    {
                      row_length = 3;
                    }
                  else
                    {
                      std::cerr << "PRISMS-PF ERROR: User-defined constant "
                                   "tensor does not have the correct number of "
                                   "elements, matrices must be 2x2 or 3x3."
                                << std::endl;
                      abort();
                    }

                  dealii::Tensor<2, dim> temp;
                  for (unsigned int i = 0; i < dim; i++)
                    {
                      for (unsigned int j = 0; j < dim; j++)
                        {
                          temp[i][j] = dealii::Utilities::string_to_double(
                            model_constants_strings.at(i * row_length + j));
                        }
                    }
                  model_constants.push_back(temp);
                }
            }
          else if (boost::iequals(model_constants_type_strings.at(1), "elastic") &&
                   boost::iequals(model_constants_type_strings.at(2), "constants"))
            {
              unsigned int num_elements = model_constants_strings.size() - 1;

              // Strip parentheses from the input, counting how many rows there
              // are
              unsigned int open_parentheses  = 0;
              unsigned int close_parentheses = 0;
              for (unsigned int element = 0; element < num_elements; element++)
                {
                  std::size_t index = 0;
                  while (index != std::string::npos)
                    {
                      index = model_constants_strings.at(element).find("(");
                      if (index != std::string::npos)
                        {
                          model_constants_strings.at(element).erase(index, 1);
                          open_parentheses++;
                        }
                    }
                  index = 0;
                  while (index != std::string::npos)
                    {
                      index = model_constants_strings.at(element).find(")");
                      if (index != std::string::npos)
                        {
                          model_constants_strings.at(element).erase(index, 1);
                          close_parentheses++;
                        }
                    }
                }
              if (open_parentheses != close_parentheses)
                {
                  std::cerr << "PRISMS-PF ERROR: User-defined elastic constant "
                               "list does not have the same number of open and "
                               "close parentheses."
                            << std::endl;
                  abort();
                }

              // Load in the elastic constants as a vector
              std::vector<double> temp_elastic_constants;
              for (unsigned int i = 0; i < num_elements; i++)
                {
                  temp_elastic_constants.push_back(
                    dealii::Utilities::string_to_double(model_constants_strings.at(i)));
                }

              std::string elastic_const_symmetry = model_constants_type_strings.at(0);
              dealii::Tensor<2, 2 *dim - 1 + dim / 3> temp =
                get_Cij_tensor(temp_elastic_constants, elastic_const_symmetry);
              model_constants.push_back(temp);
            }
          else
            {
              std::cerr << "PRISMS-PF ERROR: Only user-defined constant "
                           "tensors may have multiple elements."
                        << std::endl;
              abort();
            }
        }
    }
}

// ==========================================================================================
// Method to build the elasticity tensor from a list of elastic constants
// ==========================================================================================
template <int dim>
dealii::Tensor<2, 2 * dim - 1 + dim / 3>
userInputParameters<dim>::get_Cij_tensor(std::vector<double> elastic_constants,
                                         const std::string   elastic_const_symmetry) const
{
  // First set the material model
  elasticityModel mat_model = ISOTROPIC;
  if (elastic_const_symmetry == "isotropic")
    {
      mat_model = ISOTROPIC;
    }
  else if (elastic_const_symmetry == "transverse")
    {
      mat_model = TRANSVERSE;
    }
  else if (elastic_const_symmetry == "orthotropic")
    {
      mat_model = ORTHOTROPIC;
    }
  else if (elastic_const_symmetry == "anisotropic")
    {
      mat_model = ANISOTROPIC;
    }
  else
    {
      // Should change to an exception
      std::cerr << "Elastic material model is invalid, please use isotropic, "
                   "transverse, orthotropic, or anisotropic"
                << std::endl;
    }

  // If the material model is anisotropic for a 2D calculation but the elastic
  // constants are given for a 3D calculation, change the elastic constant
  // vector to the 2D form
  if ((mat_model == ANISOTROPIC) && (dim == 2) && elastic_constants.size() == 21)
    {
      std::vector<double> elastic_constants_temp = elastic_constants;
      elastic_constants.clear();
      std::vector<unsigned int> indices_2D = {0, 1, 5, 6, 10, 14};
      for (unsigned int i = 0; i < indices_2D.size(); i++)
        {
          elastic_constants.push_back(elastic_constants_temp.at(indices_2D.at(i)));
        }
    }

  dealii::ConditionalOStream pcout(std::cout,
                                   dealii::Utilities::MPI::this_mpi_process(
                                     MPI_COMM_WORLD) == 0);

  return getCIJMatrix(mat_model, elastic_constants, pcout);
}

// implementation of various models of anisotropy for
// st.venant-kirchoff material model of elasticity

// #include <deal.II/base/conditional_ostream.h>

// Each material model is characterized by the number of independent
// constants required to characterize its elasticity tensor in the Voigt
// notation

// 3D models:
// ISOTROPIC - 2 constants [E, nu], where E-modulus and nu-poisson's ratio
// TRANSVERSE- 5 constants [C11 C33 C44 C12 C13]
// ORTHOTROPIC- 9 constants [C11 C22 C33 C44 C55 C66 C12 C13 C23]
// ANISOTROPIC- 21 constants [C11 C22 C33 C44 C55 C66 C12 C13 C14 C15
// C16 C23 C24 C25 C26 C34 C35 C36 C45 C46 C56]

// 2D models:
// ISOTROPIC- 2 constants [E, nu] (Plane Strain)
// ANISOTROPIC- 6 constants [C11 C22 C66 C12 C16 C26] (Plane Strain)

// 1D models:
// ISOTROPIC- 1 constant [E]

template <int dim>
dealii::Tensor<2, 2 * dim - 1 + dim / 3>
userInputParameters<dim>::getCIJMatrix(const elasticityModel       model,
                                       const std::vector<double>   constants,
                                       dealii::ConditionalOStream &pcout) const
{
  // CIJ.fill(0.0);
  dealii::Tensor<2, 2 * dim - 1 + dim / 3> CIJ;

  pcout << "Reading material model:";
  switch (dim)
    {
      case 1:
        {
          pcout << " 1D ";
          // 1D models
          switch (model)
            {
              case ISOTROPIC:
                {
                  pcout << " ISOTROPIC \n";
                  CIJ[0][0] = constants[0];
                  break;
                }
              default:
                {
                  std::cout << "\nelasticityModels: Supported models in 1D - "
                               "ISOTROPIC\n";
                  std::cout << "See /src/elasticityModels.h\n";
                  exit(-1);
                }
            }
          break;
        }
      case 2:
        {
          pcout << " 2D ";
          // 2D models
          switch (model)
            {
              case ISOTROPIC:
                {
                  pcout << " ISOTROPIC \n";
                  double E = constants[0], nu = constants[1];
                  double mu     = E / (2 * (1 + nu)),
                         lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
                  CIJ[0][0]     = lambda + 2 * mu;
                  CIJ[1][1]     = lambda + 2 * mu;
                  CIJ[2][2]     = mu;
                  CIJ[0][1] = CIJ[1][0] = lambda;
                  break;
                }
              case ANISOTROPIC:
                {
                  pcout << " ANISOTROPIC \n";
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
                  std::cout << "\nelasticityModels: Supported models in 2D - "
                               "ISOTROPIC/ANISOTROPIC\n";
                  std::cout << "See /src/elasticityModels.h\n";
                  exit(-1);
                }
            }
          break;
        }
      case 3:
        {
          pcout << " 3D ";
          // 3D models
          switch (model)
            {
              case ISOTROPIC:
                {
                  pcout << " ISOTROPIC \n";
                  double E = constants[0], nu = constants[1];
                  double mu     = E / (2 * (1 + nu)),
                         lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
                  CIJ[0][0]     = lambda + 2 * mu;
                  CIJ[1][1]     = lambda + 2 * mu;
                  CIJ[2][2]     = lambda + 2 * mu;
                  CIJ[3][3]     = mu;
                  CIJ[4][4]     = mu;
                  CIJ[5][5]     = mu;
                  CIJ[0][1] = CIJ[1][0] = lambda;
                  CIJ[0][2] = CIJ[2][0] = lambda;
                  CIJ[1][2] = CIJ[2][1] = lambda;
                  break;
                }
              case TRANSVERSE:
                {
                  pcout << " TRANSVERSE \n";
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
                  pcout << " ORTHOTROPIC \n";
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
                  pcout << " ANISOTROPIC \n";
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
                  std::cout << "\nelasticityModels: Supported models in 3D - "
                               "ISOTROPIC/TRANSVERSE/ORTHOTROPIC/ANISOTROPIC\n";
                  std::cout << "See /src/elasticityModels.h\n";
                  exit(-1);
                }
            }
          break;
        }
      default:
        {
          std::cout << "\nelasticityModels: DIM is not 1/2/3\n";
          exit(-1);
        }
    }
  // print CIJ to terminal
  pcout << "Elasticity matrix (Voigt notation):\n";
  char buffer[100];
  for (unsigned int i = 0; i < 2 * dim - 1 + dim / 3; i++)
    {
      for (unsigned int j = 0; j < 2 * dim - 1 + dim / 3; j++)
        {
          snprintf(buffer, sizeof(buffer), "%8.3e ", CIJ[i][j]);
          pcout << buffer;
        }
      pcout << "\n";
    }
  pcout << "\n";
  return CIJ;
}

// Template instantiations
#include "../../include/userInputParameters_template_instantiations.h"
