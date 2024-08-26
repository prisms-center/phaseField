#include "../../include/userInputParameters.h"

// ==========================================================================================
// Method to extract the list of boundary conditions
// ==========================================================================================
template <int dim>
void
userInputParameters<dim>::load_BC_list(std::vector<std::string> list_of_BCs)
{
  // Load the BC information from the strings into a varBCs object
  // Move this to a new method and write a unit test for it!!!!

  std::vector<std::string> temp;

  for (unsigned int i = 0; i < list_of_BCs.size(); i++)
    {
      varBCs<dim> newBC;
      temp = dealii::Utilities::split_string_list(list_of_BCs[i]);

      // If there is only one BC listed, make another dim*2-1 copies of it so
      // that the same BC is applied for all boundaries
      if (temp.size() == 1)
        {
          for (unsigned int boundary = 0; boundary < (dim * 2 - 1); boundary++)
            {
              temp.push_back(temp[0]);
            }
        }

      // Load the BC for each boundary into 'newBC'
      for (unsigned int i = 0; i < (2 * dim); i++)
        {
          if (boost::iequals(temp[i], "NATURAL"))
            {
              newBC.var_BC_type.push_back(NATURAL);
              newBC.var_BC_val.push_back(0.0);
            }
          else if (boost::iequals(temp[i], "PERIODIC"))
            {
              newBC.var_BC_type.push_back(PERIODIC);
              newBC.var_BC_val.push_back(0.0);
            }
          else if (boost::iequals(temp[i], "NON_UNIFORM_DIRICHLET"))
            {
              newBC.var_BC_type.push_back(NON_UNIFORM_DIRICHLET);
              newBC.var_BC_val.push_back(0.0);
            }
          else if (boost::iequals(temp[i].substr(0, 9), "DIRICHLET"))
            {
              newBC.var_BC_type.push_back(DIRICHLET);
              std::string dirichlet_val = temp[i].substr(10, temp[i].size());
              dirichlet_val             = dealii::Utilities::trim(dirichlet_val);
              newBC.var_BC_val.push_back(
                dealii::Utilities::string_to_double(dirichlet_val));
            }
          else if (boost::iequals(temp[i].substr(0, 7), "NEUMANN"))
            {
              newBC.var_BC_type.push_back(NEUMANN);
              std::string neumann_val = temp[i].substr(8, temp[i].size());
              neumann_val             = dealii::Utilities::trim(neumann_val);
              newBC.var_BC_val.push_back(
                dealii::Utilities::string_to_double(neumann_val));
            }
          else
            {
              std::cout << temp[i].substr(0, 8) << std::endl;
              std::cout << "Error: Boundary conditions specified improperly."
                        << std::endl;
              abort();
            }
        }
      BC_list.push_back(newBC);

      // Validate input using something like this:
      // try{
      //     if ((BC_type_dim1_min == "PERIODIC") && (BC_type_dim1_max !=
      //     "PERIODIC")){
      //         throw 0;
      //     }
      //     if ((BC_type_dim2_min == "PERIODIC") && (BC_type_dim2_max !=
      //     "PERIODIC")){
      //         throw 0;
      //     }
      //     if ((BC_type_dim3_min == "PERIODIC") && (BC_type_dim3_max !=
      //     "PERIODIC")){
      //         throw 0;
      //     }
      // }
      // catch (int e){
      //     if (e == 0){
      //         std::cout << "Error: For periodic BCs, both faces for a given
      //         direction must be set as periodic. "
      //                 "Please check the BCs that are set in ICs_and_BCs.h."
      //                 << std::endl;
      //     }
      //     abort();
      // }
    }
}

// Template instantiations
#include "../../include/userInputParameters_template_instantiations.h"
