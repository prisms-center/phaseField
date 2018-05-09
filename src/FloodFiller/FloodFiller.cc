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

    // The flood fill loop
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

        ++di;
    }

    if (grain_sets.back().getVertexList().size() == 0){
        grain_sets.pop_back();
    }

    // Get a global list of grains from the various local lists
    communicateGrainSets(grain_sets);

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

                std::vector<dealii::Point<dim>> vertex_list;
                for (unsigned int v=0; v< dealii::Utilities::fixed_power<dim>(2.0); v++){
                    vertex_list.push_back(di->vertex(v));
                }
                grain_sets.back().addVertexList(vertex_list);

                // Call recursiveFloodFill on the element's neighbors
                for (unsigned int n=0; n<2*dim; n++){
                    recursiveFloodFill<typename DoFHandler<dim>::active_cell_iterator>(di->neighbor(n), di_end, solution_field, threshold,  grain_index, grain_sets, grain_assigned);
                }
            }
        }
    }
}

// =================================================================================
// Generate global list of the grains, merging grains split between multiple processors
// =================================================================================
template <int dim, int degree>
void FloodFiller<dim, degree>::communicateGrainSets(std::vector<GrainSet<dim>> & grain_sets){

    //MPI INITIALIZATON
	int numProcs=dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
	int thisProc=dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

	if (numProcs > 1) {
		// Cycle through each processor, sending and receiving, to append the list of grains
		for (int proc_index=0; proc_index < numProcs-1; proc_index++){
			if (thisProc == proc_index){
				sendUpdate(thisProc+1, grain_sets);
			}
			else if (thisProc == proc_index+1){
				receiveUpdate(thisProc-1, grain_sets);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		// The final processor now has all of the grains
		// Check for grains that are split across processors
		if (thisProc == numProcs-1){
			mergeSplitGrains(grain_sets);
		}

		// The final processor now has the final list of the grains, broadcast it to all the other processors
		broadcastUpdate(numProcs-1, thisProc, grain_sets);
	}

}

// =================================================================================
// Sends the list of new nuclei to the next processor
// =================================================================================
template <int dim, int degree>
void FloodFiller<dim, degree>::sendUpdate (int procno, std::vector<GrainSet<dim>> & grain_sets) const
{
    int num_grains=grain_sets.size();
    //MPI SECTION TO SEND INFORMATION TO THE PROCESSOR procno
    //Sending local no. of grains
    MPI_Send(&num_grains, 1, MPI_INT, procno, 0, MPI_COMM_WORLD);
    if (num_grains > 0){

        std::vector<unsigned int> s_num_elements;
        std::vector<double> s_vertices;

        for (unsigned int g=0; g<grain_sets.size(); g++){

            std::vector<std::vector<dealii::Point<dim>>> vertex_list = grain_sets[g].getVertexList();
            s_num_elements.push_back(vertex_list.size());

            for (unsigned int c=0; c<s_num_elements[g]; c++){
                for (unsigned int v=0; v<dealii::Utilities::fixed_power<dim>(2.0); v++){
                    for (unsigned int d=0; d<dim; d++){
                        s_vertices.push_back(vertex_list[c][v][d]);
                    }
                }
            }
        }

        unsigned int num_vertices = 0;
        for (unsigned int g=0; g<grain_sets.size(); g++){
            num_vertices += s_num_elements[g] * dealii::Utilities::fixed_power<dim>(2.0) * (double)dim;
        }

        MPI_Send(&s_num_elements[0], num_grains, MPI_UNSIGNED, procno, 1, MPI_COMM_WORLD);
        MPI_Send(&s_vertices[0], num_vertices, MPI_DOUBLE, procno, 2, MPI_COMM_WORLD);
    }
}

// =================================================================================
// Recieves the list of new nuclei from the previous processor
// =================================================================================
template <int dim, int degree>
void FloodFiller<dim, degree>::receiveUpdate (int procno, std::vector<GrainSet<dim>> & grain_sets) const
{
    //MPI PROCEDURE TO RECIEVE INFORMATION FROM ANOTHER PROCESSOR AND UPDATE LOCAL NUCLEI INFORMATION
    int num_grains = 0;
    MPI_Recv(&num_grains, 1, MPI_INT, procno, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (num_grains > 0){

        std::vector<unsigned int> r_num_elements(num_grains, 0);
        MPI_Recv(&r_num_elements[0], num_grains, MPI_UNSIGNED, procno, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        unsigned int num_vertices = 0;
        for (int g=0; g<num_grains; g++){
            num_vertices += r_num_elements[g] * dealii::Utilities::fixed_power<dim>(2.0) * (double)dim;
        }

        std::vector<double> r_vertices(num_vertices, 0.0);
        MPI_Recv(&r_vertices[0], num_vertices, MPI_DOUBLE, procno, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Now add these new grains to the local grain list
        for (int g=0; g<num_grains; g++){
            GrainSet<dim> new_grain_set;
            new_grain_set.setGrainIndex(grain_sets.size()+1);
            for (unsigned int c=0; c<r_num_elements.at(g); c++){
                std::vector<dealii::Point<dim>> verts;
                for (unsigned int v=0; v < dealii::Utilities::fixed_power<dim>(2.0); v++){
                    double coords[dim];
                    for (unsigned int d=0; d<dim; d++){
                        coords[d] = r_vertices.front();
                        r_vertices.erase(r_vertices.begin());
                    }
                    dealii::Tensor<1,dim> tensor_coords(coords);
                    dealii::Point<dim> vert(tensor_coords);
                    verts.push_back(vert);
                }
                new_grain_set.addVertexList(verts);
            }
            grain_sets.push_back(new_grain_set);
        }

    }
}

// =================================================================================
// Broadcast the final list of new nuclei from the last processor to the rest
// =================================================================================
template <int dim, int degree>
void FloodFiller<dim, degree>::broadcastUpdate (int broadcastProc, int thisProc, std::vector<GrainSet<dim>> & grain_sets) const
{
    //MPI PROCEDURE TO SEND THE LIST OF GRAINS FROM ONE PROCESSOR TO ALL THE OTHERS
    int num_grains=grain_sets.size();
    MPI_Bcast(&num_grains, 1, MPI_INT, broadcastProc, MPI_COMM_WORLD);
    if (num_grains > 0){

        //Creating vectors of each quantity in nuclei. Each numbered acording to the tags used for MPI_Send/MPI_Recv
    	unsigned int initial_vec_size;
    	if (thisProc == broadcastProc){
    		initial_vec_size = 0;
    	}
    	else{
    		initial_vec_size = num_grains;
    	}

        std::vector<unsigned int> num_elements(initial_vec_size, 0);
        std::vector<double> vertices;

        if (thisProc == broadcastProc){
            for (unsigned int g=0; g<grain_sets.size(); g++){

                std::vector<std::vector<dealii::Point<dim>>> vertex_list = grain_sets[g].getVertexList();
                num_elements.push_back(vertex_list.size());

                for (unsigned int c=0; c<num_elements[g]; c++){
                    for (unsigned int v=0; v<dealii::Utilities::fixed_power<dim>(2.0); v++){
                        for (unsigned int d=0; d<dim; d++){
                            vertices.push_back(vertex_list[c][v][d]);
                        }
                    }
                }
            }
        }

        MPI_Bcast(&num_elements[0], num_grains, MPI_UNSIGNED, broadcastProc, MPI_COMM_WORLD);

        unsigned int num_vertices = 0;
        for (int g=0; g<num_grains; g++){
            num_vertices += num_elements[g] * dealii::Utilities::fixed_power<dim>(2.0) * (double)dim;
        }

        if (thisProc != broadcastProc){
            vertices.assign(num_vertices,0.0);
        }

        MPI_Bcast(&vertices[0], num_vertices, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);

        // Now add these new grains to the local grain list
        if (thisProc != broadcastProc){
            grain_sets.clear();
            for (int g=0; g<num_grains; g++){
                GrainSet<dim> new_grain_set;
                new_grain_set.setGrainIndex(grain_sets.size()+1);
                for (unsigned int c=0; c<num_elements.at(g); c++){
                    std::vector<dealii::Point<dim>> verts;
                    for (unsigned int v=0; v < dealii::Utilities::fixed_power<dim>(2.0); v++){
                        double coords[dim];
                        for (unsigned int d=0; d<dim; d++){
                            coords[d] = vertices.front();
                            vertices.erase(vertices.begin());
                        }
                        dealii::Tensor<1,dim> tensor_coords(coords);
                        dealii::Point<dim> vert(tensor_coords);
                        verts.push_back(vert);
                    }
                    new_grain_set.addVertexList(verts);
                }
                grain_sets.push_back(new_grain_set);
            }
        }
    }
}

template <int dim, int degree>
void FloodFiller<dim, degree>::mergeSplitGrains (std::vector<GrainSet<dim>> & grain_sets) const
{
    // Loop though each vertex in the base grain "g"
    for (unsigned int g=0; g<grain_sets.size(); g++){

        std::vector<std::vector<dealii::Point<dim>>> vertex_list = grain_sets[g].getVertexList();

        for (unsigned int c = 0; c < vertex_list.size(); c++){
            for (unsigned int v = 0; v < dealii::Utilities::fixed_power<dim>(2.0); v++){

                // Now cycle through the other grains to find overlapping elements
                for (unsigned int g_other=g+1; g_other<grain_sets.size(); g_other++){
                    bool matching_vert = false;

                    std::vector<std::vector<dealii::Point<dim>>> vertex_list_other = grain_sets[g_other].getVertexList();

                    for (unsigned int c_other = 0; c_other < vertex_list_other.size(); c_other++){
                        for (unsigned int v_other = 0; v_other < dealii::Utilities::fixed_power<dim>(2.0); v_other++){
                            // Check if the vertices match
                            if (vertex_list[c][v] == vertex_list_other[c_other][v_other]){
                                matching_vert = true;
                                break;
                            }
                            if (matching_vert){
                                break;
                            }
                        }
                        if (matching_vert){
                            break;
                        }
                    }

                    // If this grain has a vertex shared with the base grain, move the vertices to the base grain and delete the other grain
                    if (matching_vert){
                        for (unsigned int c_other = 0; c_other < vertex_list_other.size(); c_other++){
                            grain_sets[g].addVertexList(vertex_list_other.at(c_other));
                        }
                        grain_sets.erase(grain_sets.begin()+g_other);
                    }
                }
            }
        }
    }
}
