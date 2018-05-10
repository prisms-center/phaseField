#include "../../include/OrderParameterRemapper.h"

template <int dim>
void OrderParameterRemapper<dim>::remap(
    std::vector<SimplifiedGrainRepresentation<dim>> & grain_representations,
    std::vector<vectorType*> solution_fields, dealii::DoFHandler<dim> & dof_handler)
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
                    bool in_grain = false;
                    for (unsigned int v=0; v< dealii::Utilities::fixed_power<dim>(2.0); v++){
                        if (di->vertex(v).distance(grain_representations.at(g).getCenter())){
                            in_grain = true;
                            break;
                        }
                    }

                    // If it is, move the values from the old order parameter to the new order parameter and set the old order parameter to zero
                    if (in_grain){
                        std::vector<types::global_dof_index> dof_indices;
                        di->get_dof_indices(dof_indices);

                        for (unsigned int i=0; i < dof_indices.size(); i++){
                            solution_fields.at(op_new)[i] = solution_fields.at(op_old)[i];
                            solution_fields.at(op_old)[i] = 0.0;
                        }
                    }
                }
                ++di;
            }
        }
    }
}
