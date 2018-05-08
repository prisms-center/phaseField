#include "../../include/FloodFiller.h"

template <int dim, int degree>
void FloodFiller<dim, degree>::calcGrainSets(FESystem<dim> & fe, dealii::DoFHandler<dim> &dof_handler, vectorType* solution_field, double threshold, std::vector<GrainSet<dim>> & grain_sets){

    unsigned int grain_index = 0;

    // Loop through the whole mesh and set the user flags to false (so everything is considered unmarked)
    typename DoFHandler<dim>::active_cell_iterator di = dof_handler.begin_active();
    while (di != dof_handler.end())
    {
        di->clear_user_flag();
        ++di;
    }

    GrainSet<dim> grain_set;
    grain_sets.push_back(grain_set);
    grain_sets.back().setGrainIndex(grain_index);

    di = dof_handler.begin_active();
    while (di != dof_handler.end())
    {
        bool grain_assigned = false;
        recursiveFloodFill<typename DoFHandler<dim>::active_cell_iterator>(di, dof_handler.end(), solution_field, threshold,  grain_index, grain_sets, grain_assigned);


        if (grain_assigned){
            grain_index++;
            GrainSet<dim> grain_set;
            grain_sets.push_back(grain_set);
            grain_sets.back().setGrainIndex(grain_index);
        }

        for (unsigned int g=0; g<grain_sets.size(); g++ ){
            std::cout << "Grain index: " << grain_sets.at(g).getGrainIndex() << " Vertex list length: " << grain_sets.at(g).getVertexList().size() << std::endl;
        }


        ++di;
    }

}

template <int dim, int degree>
template <typename T>
void FloodFiller<dim, degree>::recursiveFloodFill(T di, T di_end, vectorType* solution_field, double threshold, unsigned int & grain_index, std::vector<GrainSet<dim>> & grain_sets, bool & grain_assigned){

    if (di != di_end and di->is_locally_owned()){

        // Check if the cell has been marked yet
        bool cellMarked = di->user_flag_set();

        if (!cellMarked){

            di->set_user_flag();

            FEValues<dim> fe_values (*fe, quadrature, update_values);
            std::vector<double> var_values(num_quad_points);
            std::vector<dealii::Point<dim> > q_point_list(num_quad_points);

            // Get the average value for the element
            fe_values.reinit(di);
            fe_values.get_function_values(*solution_field, var_values);

            double ele_val = 0.0;
            for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                for (unsigned int i=0; i<dofs_per_cell; ++i){
                    ele_val += fe_values.shape_value (i, q_point)*var_values[q_point]*quadrature.weight(q_point);
                }
            }


            if (ele_val > threshold){
                grain_assigned = true;

                //grain_sets.back().setGrainIndex(grain_index);
                std::vector<dealii::Point<dim>> vertex_list;
                for (unsigned int v=0; v< dealii::Utilities::fixed_power<dim>(2.0); v++){
                    vertex_list.push_back(di->vertex(v));
                }
                grain_sets.back().addVertexList(vertex_list);

                std::cout << "(in) grain index: " << grain_sets.back().getGrainIndex() << " Vertex list length: " << grain_sets.back().getVertexList().size() << std::endl;

                for (unsigned int n=0; n<2*dim; n++){
                    recursiveFloodFill<typename DoFHandler<dim>::active_cell_iterator>(di->neighbor(n), di_end, solution_field, threshold,  grain_index, grain_sets, grain_assigned);
                }
            }
        }
    }
}
