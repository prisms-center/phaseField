#include "../../include/OrderParameterRemapper.h"

template <int dim>
void OrderParameterRemapper<dim>::remap(
    std::vector<SimplifiedGrainRepresentation<dim>> & grain_representations,
    std::vector<vectorType*> & solution_fields, dealii::DoFHandler<dim> & dof_handler, unsigned int dofs_per_cell)
{
    for (unsigned int g=0; g < grain_representations.size(); g++){
        if (grain_representations.at(g).getOrderParameterId() != grain_representations.at(g).getOldOrderParameterId()){

            typename DoFHandler<dim>::active_cell_iterator di = dof_handler.begin_active();

            while (di != dof_handler.end())
            {
                if (di->is_locally_owned()){
                    unsigned int op_new = grain_representations.at(g).getOrderParameterId();
                    unsigned int op_old = grain_representations.at(g).getOldOrderParameterId();

                    // Check if the cell is within the simplified grain representation
                    bool in_grain = true;
                    for (unsigned int v=0; v< dealii::GeometryInfo<dim>::vertices_per_cell; v++){
                        if (di->vertex(v).distance(grain_representations.at(g).getCenter()) > grain_representations.at(g).getRadius()){
                            in_grain = false;
                            break;
                        }
                    }

                    // If it is, move the values from the old order parameter to the new order parameter and set the old order parameter to zero
                    if (in_grain){
                        std::vector<types::global_dof_index> dof_indices(dofs_per_cell,0);
                        di->get_dof_indices(dof_indices);

                        for (unsigned int i=0; i < dof_indices.size(); i++){
                            //std::cout << "b" << solution_fields.at(op_old)->local_element(dof_indices.at(i)) << " " << solution_fields.at(op_new)->local_element(dof_indices.at(i)) <<  std::endl;

                            // Something isn't right in these three lines, the problem is that I'm modifying the fields in place so sometimes I'm not getting the correct values when I assign the value from the new OP to the old OP. I need to break this up into stages or create a temp vector.
                            solution_fields.at(op_new)->local_element(dof_indices.at(i)) = solution_fields.at(op_old)->local_element(dof_indices.at(i));

                            //solution_fields.at(op_old)->local_element(dof_indices.at(i)) = 0.0;

                            //std::cout << solution_fields.at(op_old)->local_element(dof_indices.at(i)) << "\t";

                            //std::cout << "a" << solution_fields.at(op_old)->local_element(dof_indices.at(i)) << " " << solution_fields.at(op_new)->local_element(dof_indices.at(i)) <<  std::endl;
                        }
                        //std::cout << "\n";
                    }
                }
                ++di;
            }

            di = dof_handler.begin_active();

            while (di != dof_handler.end())
            {
                if (di->is_locally_owned()){
                    unsigned int op_new = grain_representations.at(g).getOrderParameterId();
                    unsigned int op_old = grain_representations.at(g).getOldOrderParameterId();

                    // Check if the cell is within the simplified grain representation
                    bool in_grain = true;
                    for (unsigned int v=0; v< dealii::GeometryInfo<dim>::vertices_per_cell; v++){
                        if (di->vertex(v).distance(grain_representations.at(g).getCenter()) > grain_representations.at(g).getRadius()){
                            in_grain = false;
                            break;
                        }
                    }

                    // If it is, set the old order parameter to zero
                    if (in_grain){
                        std::vector<types::global_dof_index> dof_indices(dofs_per_cell,0);
                        di->get_dof_indices(dof_indices);

                        for (unsigned int i=0; i < dof_indices.size(); i++){

                            solution_fields.at(op_old)->local_element(dof_indices.at(i)) = 0.0;

                        }
                    }
                }
                ++di;
            }
        }
    }
}
