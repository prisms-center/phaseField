// ------------------------------------------------------------------------
// Method to extract the list of boundary conditions from input parameters
// and store them in BC_list object
// ------------------------------------------------------------------------

#include <deal.II/base/exceptions.h>

#include "../../include/userInputParameters.h"

template <int dim>
void
userInputParameters<dim>::load_BC_list(const std::vector<std::string> &list_of_BCs)
{
  // Loop over the list of boundary conditions specified in parameters
  // and provided in the input list_of_BCs. Process the BCs and place
  // them into the vector BC_list
  std::vector<std::string> temp;
  for (const auto &boundary_condition : list_of_BCs)
    {
      // Ensure all variables have BCs specified in parameters.prm
      AssertThrow(!boundary_condition.empty(),
                  dealii::ExcMessage(std::string("Boundary condition not specified.")));

      varBCs<dim> newBC;
      temp = dealii::Utilities::split_string_list(boundary_condition);

      // If there is only one BC listed, make another 2*dim-1 copies of it so that
      // the same BC is applied for all boundaries
      if (temp.size() == 1)
        {
          temp.resize(2 * dim, temp[0]);
        }

      // Load the BC for each boundary into 'newBC'.
      for (unsigned int j = 0; j < (2 * dim); j++)
        {
          if (boost::iequals(temp[j], "NATURAL"))
            {
              newBC.var_BC_type.push_back(NATURAL);
              newBC.var_BC_val.push_back(0.0);
            }
          else if (boost::iequals(temp[j], "PERIODIC"))
            {
              newBC.var_BC_type.push_back(PERIODIC);
              newBC.var_BC_val.push_back(0.0);
            }
          else if (boost::iequals(temp[j], "NON_UNIFORM_DIRICHLET"))
            {
              newBC.var_BC_type.push_back(NON_UNIFORM_DIRICHLET);
              newBC.var_BC_val.push_back(0.0);
            }
          else if (boost::iequals(temp[j].substr(0, 9), "DIRICHLET"))
            {
              newBC.var_BC_type.push_back(DIRICHLET);
              std::string dirichlet_val = temp[j].substr(10, temp[j].size());
              dirichlet_val             = dealii::Utilities::trim(dirichlet_val);
              newBC.var_BC_val.push_back(
                dealii::Utilities::string_to_double(dirichlet_val));
            }
          else if (boost::iequals(temp[j].substr(0, 7), "NEUMANN"))
            {
              newBC.var_BC_type.push_back(NEUMANN);
              std::string neumann_val = temp[j].substr(8, temp[j].size());
              neumann_val             = dealii::Utilities::trim(neumann_val);
              newBC.var_BC_val.push_back(
                dealii::Utilities::string_to_double(neumann_val));
            }
          else
            {
              std::cout << temp[j].substr(0, 8) << std::endl;
              std::cout << "Error: Boundary conditions specified improperly."
                        << std::endl;
              abort();
            }

          // If periodic BCs are used, ensure they are applied on both sides of
          // domain
          if (j % 2 == 0)
            {
              AssertThrow(!((boost::iequals(temp[j], "PERIODIC") &&
                             !boost::iequals(temp[j + 1], "PERIODIC")) ||
                            (!boost::iequals(temp[j], "PERIODIC") &&
                             boost::iequals(temp[j + 1], "PERIODIC"))),
                          dealii::ExcMessage(
                            std::string("Periodic boundary condition must be "
                                        "specified on both sides of domain")));
            }
        }
      // Append BCs for current field to total list
      BC_list.push_back(newBC);
    }
}

#include "../../include/userInputParameters_template_instantiations.h"