// Methods in MatrixFreePDE to update the list of nuclei
#include "../../include/matrixFreePDE.h"
#include "../../include/parallelNucleationList.h"
#include <random>
#include <time.h>

// =======================================================================================================
// Function called in solve to update the global list of nuclei
// =======================================================================================================
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::updateNucleiList() {

    if (userInputs.nucleation_occurs){


        if (currentIncrement % userInputs.steps_between_nucleation_attempts == 0 || currentIncrement == 1){
            computing_timer.enter_section("matrixFreePDE: nucleation");

            // Apply constraints
            for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
                constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                solutionSet[fieldIndex]->update_ghost_values();
            }

            std::vector<nucleus<dim> > new_nuclei;

            if (currentIncrement == 1){
                while (new_nuclei.size() == 0){
                    currentTime+=userInputs.dtValue*(double)userInputs.steps_between_nucleation_attempts;
                    currentIncrement+=userInputs.steps_between_nucleation_attempts;

                    while (userInputs.outputTimeStepList.size() > 0 && userInputs.outputTimeStepList[currentOutput] < currentIncrement){
                        currentOutput++;
                    }

                    while (userInputs.checkpointTimeStepList.size() > 0 && userInputs.checkpointTimeStepList[currentCheckpoint] < currentIncrement){
                        currentCheckpoint++;
                    }

                    new_nuclei = getNewNuclei();
                }
            }
            else {
                new_nuclei = getNewNuclei();
            }

            nuclei.insert(nuclei.end(),new_nuclei.begin(),new_nuclei.end());

            if (new_nuclei.size() > 0 && userInputs.h_adaptivity == true){
                refineMeshNearNuclei(new_nuclei);
            }

            computing_timer.exit_section("matrixFreePDE: nucleation");
        }
    }

}

// =======================================================================================================
// Core method to perform a nucleation check
// =======================================================================================================
template <int dim, int degree>
std::vector<nucleus<dim> > MatrixFreePDE<dim,degree>::getNewNuclei(){

    // Declare a vector of all the NEW nuclei seeded in this time step
    std::vector<nucleus<dim> > newnuclei;

    // Get list of prospective new nuclei for the local processor
    pcout << "Nucleation attempt for increment " << currentIncrement << std::endl;

    getLocalNucleiList(newnuclei);
    pcout << "nucleation attempt! " << currentTime << " " << currentIncrement << std::endl;

    // Generate global list of new nuclei and resolve conflicts between new nuclei
    parallelNucleationList<dim> new_nuclei_parallel(newnuclei);
    newnuclei = new_nuclei_parallel.buildGlobalNucleiList(userInputs.min_distance_between_nuclei, nuclei.size());

    // Final check to resolve overlap conflicts with existing precipitates
    std::vector<unsigned int> conflict_ids;
    safetyCheckNewNuclei(newnuclei, conflict_ids);

    newnuclei = new_nuclei_parallel.removeSubsetOfNuclei(conflict_ids, nuclei.size());

    return newnuclei;
}

// =================================================================================
// Get list of prospective new nuclei for the local processor
// =================================================================================
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::getLocalNucleiList(std::vector<nucleus<dim> > &newnuclei) const
{
    // Nickname for current time and time step
    double t=currentTime;
    unsigned int inc=currentIncrement;

    //QGauss<dim>  quadrature(degree+1);
    QGaussLobatto<dim>  quadrature(degree+1);
    FEValues<dim> fe_values (*(FESet[0]), quadrature, update_values|update_quadrature_points|update_JxW_values);
    const unsigned int   num_quad_points = quadrature.size();
    std::vector<std::vector<double> > var_values(userInputs.nucleation_need_value.size(),std::vector<double>(num_quad_points));
    std::vector<dealii::Point<dim> > q_point_list(num_quad_points);

    std::vector<dealii::Point<dim> > q_point_list_overlap(num_quad_points);

    typename DoFHandler<dim>::active_cell_iterator   di = dofHandlersSet_nonconst[0]->begin_active();

    // What used to be in nuc_attempt
    double rand_val;
    //Better random no. generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distr(0.0,1.0);

    //Element cycle
    while (di != this->dofHandlersSet_nonconst[0]->end())
    {
        if (di->is_locally_owned()){
            //Obtaining average element concentration by averaging over element's quadrature points
            fe_values.reinit(di);
            for (unsigned int var = 0; var < userInputs.nucleation_need_value.size(); var++){
                fe_values.get_function_values(*(solutionSet[userInputs.nucleation_need_value[var]]), var_values[var]);
            }
            q_point_list = fe_values.get_quadrature_points();

            // ---------------------------

            double element_volume = 0.0;
            dealii::Point<dim> ele_center;
            // Loop over the quadrature points to find the element volume (or area in 2D) and the average q point location
            for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                element_volume = element_volume + fe_values.JxW(q_point);
                for (unsigned int i=0; i<dim; i++)
                ele_center[i]=ele_center[i]+q_point_list[q_point](i)/((double)num_quad_points);
            }

            // Loop over each variable and each quadrature point to get the average variable value for the element
            variableValueContainer variable_values;
            for (unsigned int var = 0; var < userInputs.nucleation_need_value.size(); var++){
                double ele_val = 0.0;
                for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                    ele_val = ele_val + var_values[var][q_point]*fe_values.JxW(q_point);
                }
                ele_val /= element_volume;
                variable_values.set(userInputs.nucleation_need_value[var],ele_val);
            }

            // Loop through each nucleating order parameter
            for (unsigned int i = 0; i < userInputs.nucleating_variable_indices.size(); i++){
                unsigned int variable_index = userInputs.nucleating_variable_indices.at(i);

                //Compute random no. between 0 and 1 (new method)
                rand_val=distr(gen);
                //Nucleation probability
                double Prob=getNucleationProbability(variable_values,element_volume,ele_center,variable_index);

                // ----------------------------

                if (rand_val <= Prob){

                    //Initializing random vector in "dim" dimensions
                    std::vector<double> randvec(dim,0.0);
                    dealii::Point<dim> nuc_ele_pos;

                    //Finding coordinates of quadrature point closest to and furthest away from the origin
                    std::vector<double> ele_origin(dim);
                    for (unsigned int i=0; i<dim; i++)
                    ele_origin[i] = q_point_list[0](i);
                    std::vector<double> ele_max(dim);
                    for (unsigned int i=0; i<dim; i++)
                    ele_max[i] = q_point_list[0](i);
                    for (unsigned int i=0; i<dim; i++){
                        for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                            for (unsigned int i=0; i<dim; i++){
                                if (q_point_list[q_point](i) < ele_origin[i])
                                ele_origin[i]=q_point_list[q_point](i);
                                if (q_point_list[q_point](i) > ele_max[i])
                                ele_max[i]=q_point_list[q_point](i);
                            }
                        }
                    }

                    //Find a random point within the element
                    for (unsigned int i=0; i<dim; i++){
                        randvec[i]=distr(gen);
                        nuc_ele_pos[i]=ele_origin[i] + (ele_max[i]-ele_origin[i])*randvec[i];
                    }

                    //Make sure point is in safety zone
                    bool insafetyzone = true;
                    for (unsigned int j=0; j < dim; j++){
                        bool periodic_j = (userInputs.BC_list[1].var_BC_type[2*j]==PERIODIC);
                        bool insafetyzone_j = (periodic_j || ((nuc_ele_pos[j] > userInputs.nucleation_parameters_list[i].no_nucleation_border_thickness) && (nuc_ele_pos[j] < userInputs.domain_size[j]-userInputs.nucleation_parameters_list[i].no_nucleation_border_thickness)));
                        insafetyzone = insafetyzone && insafetyzone_j;
                    }

                    if (insafetyzone){

                        // Check to see if the order parameter anywhere within the element is above the threshold
                        bool anyqp_OK = false;
                        for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                            double sum_op = 0.0;
                            for (unsigned int var = 0; var < userInputs.nucleation_need_value.size(); var++){
                                for (unsigned int op = 0; op < userInputs.nucleating_variable_indices.size(); op++){
                                    if (userInputs.nucleation_need_value[var] == userInputs.nucleating_variable_indices[op]){
                                        sum_op += var_values[var][q_point];
                                    }
                                }
                            }
                            if (sum_op < userInputs.nucleation_order_parameter_cutoff){
                                anyqp_OK =true;
                            }
                        }

                        if (anyqp_OK){
                            // Pick the order parameter (not needed anymore since the probability is now calculated on a per OP basis)
                            /*
                            std::random_device rd2;
                            std::mt19937 gen2(rd2());
                            std::uniform_int_distribution<unsigned int> int_distr(0,userInputs.nucleating_variable_indices.size()-1);
                            unsigned int op_for_nucleus = userInputs.nucleating_variable_indices[int_distr(gen2)];
                            std::cout << "Nucleation order parameter: " << op_for_nucleus << " " << rand_val << std::endl;
                            */

                            //Add nucleus to prospective list
                            std::cout << "Prospective nucleation event. Nucleus no. " << nuclei.size()+1 << std::endl;
                            std::cout << "Nucleus center: " << nuc_ele_pos << std::endl;
                            std::cout << "Nucleus order parameter: " << variable_index << std::endl;
                            nucleus<dim>* temp = new nucleus<dim>;
                            temp->index=nuclei.size();
                            temp->center=nuc_ele_pos;
                            temp->semiaxes = userInputs.nucleation_parameters_list[i].semiaxes;
                            temp->seededTime=t;
                            temp->seedingTime = userInputs.nucleation_parameters_list[i].hold_time;
                            temp->seedingTimestep = inc;
                            temp->orderParameterIndex = variable_index;
                            newnuclei.push_back(*temp);
                        }
                    }
                }
            }
        }
        // Increment the cell iterators
        ++di;
    }
}

// =======================================================================================================
// Making sure all new nuclei from complete prospective list do not overlap with existing precipitates
// =======================================================================================================
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::safetyCheckNewNuclei(std::vector<nucleus<dim> > newnuclei, std::vector<unsigned int> &conflict_ids)
{
    //QGauss<dim>  quadrature(degree+1);
    QGaussLobatto<dim>  quadrature(degree+1);
    FEValues<dim> fe_values (*(FESet[0]), quadrature, update_values|update_quadrature_points|update_JxW_values);
    const unsigned int   num_quad_points = quadrature.size();
    std::vector<std::vector<double> > op_values(userInputs.nucleating_variable_indices.size(),std::vector<double>(num_quad_points));
    std::vector<dealii::Point<dim> > q_point_list(num_quad_points);

    //Nucleus cycle
    for (typename std::vector<nucleus<dim> >::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){
        bool isClose=false;

        unsigned int nucleation_parameters_list_index;
        for (unsigned int j=0; j<userInputs.nucleation_parameters_list.size(); j++ ){
            if (userInputs.nucleation_parameters_list[j].var_index == thisNuclei->orderParameterIndex){
                nucleation_parameters_list_index = j;
            }
        }

        //Element cycle
	    typename DoFHandler<dim>::active_cell_iterator   di = dofHandlersSet_nonconst[0]->begin_active();
        while (di != dofHandlersSet_nonconst[0]->end())
        {
            if (di->is_locally_owned()){
                fe_values.reinit(di);
                for (unsigned int var = 0; var < userInputs.nucleating_variable_indices.size(); var++){
                	fe_values.get_function_values(*(solutionSet[userInputs.nucleating_variable_indices[var]]), op_values[var]);
                }
                q_point_list = fe_values.get_quadrature_points();

                //Quadrature points cycle
                for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                    // Calculate the ellipsoidal distance to the center of the nucleus
                    double weighted_dist = 0.0;
                    for (unsigned int i=0; i<dim; i++){
                        double shortest_edist = thisNuclei->center(i) - q_point_list[q_point](i);
                        bool periodic_i = (userInputs.BC_list[1].var_BC_type[2*i]==PERIODIC);
                        if (periodic_i){
                            double domsize = userInputs.domain_size[i];
                            shortest_edist = shortest_edist-round(shortest_edist/domsize)*domsize;
                        }

                        double temp = shortest_edist/(userInputs.nucleation_parameters_list[nucleation_parameters_list_index].freeze_semiaxes[i]);
                        weighted_dist += temp*temp;
                    }
                    if (weighted_dist < 1.0){
                    	double sum_op = 0.0;
                    	for (unsigned int num_op = 0; num_op < userInputs.nucleating_variable_indices.size(); num_op++){
                    		sum_op += op_values[num_op][q_point];
                    	}
                        if (sum_op > 0.1){
                            isClose=true;
                            std::cout << "Attempted nucleation failed due to overlap w/ existing particle!!!!!!"  << std::endl;
                            conflict_ids.push_back(thisNuclei->index);
                            break;
                        }
                    }
                }
                if (isClose) break;
            }
        // Increment the cell iterators
        ++di;
        }
    }
}

// =================================================================================
// Refine mesh near the new nuclei
// =================================================================================
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::refineMeshNearNuclei(std::vector<nucleus<dim> > newnuclei)
{
	//QGauss<dim>  quadrature(degree+1);
	QGaussLobatto<dim>  quadrature(degree+1);
	FEValues<dim> fe_values (*(FESet[0]), quadrature, update_values|update_quadrature_points|update_JxW_values);
	const unsigned int   num_quad_points = quadrature.size();
	std::vector<dealii::Point<dim> > q_point_list(num_quad_points);

	typename Triangulation<dim>::active_cell_iterator ti;
	typename DoFHandler<dim>::active_cell_iterator   di;

	unsigned int numDoF_preremesh = totalDOFs;

	for (unsigned int remesh_index=0; remesh_index < (userInputs.max_refinement_level-userInputs.min_refinement_level); remesh_index++){
		ti  = triangulation.begin_active();
		di = dofHandlersSet_nonconst[0]->begin_active();
		while (di != dofHandlersSet_nonconst[0]->end()){
			if (di->is_locally_owned()){

				bool mark_refine = false;

				fe_values.reinit (di);
				q_point_list = fe_values.get_quadrature_points();

				// Calculate the distance from the corner of the cell to the middle of the cell
				double diag_dist = 0.0;
				for (unsigned int i=0; i<dim; i++){
					diag_dist += (userInputs.domain_size[i]*userInputs.domain_size[i])/(userInputs.subdivisions[i]*userInputs.subdivisions[i]);
				}
				diag_dist = sqrt(diag_dist);
				diag_dist /= 2.0*pow(2.0,ti->level());

				for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
					for (typename std::vector<nucleus<dim> >::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){

                        // Calculate the ellipsoidal distance to the center of the nucleus
						double weighted_dist = 0.0;
                        dealii::Tensor<1,dim,double> shortest_edist_tensor = thisNuclei->center - q_point_list[q_point];
                        for (unsigned int i=0; i<dim; i++){
                            if (userInputs.BC_list[thisNuclei->orderParameterIndex].var_BC_type[2*i]==PERIODIC){
                                shortest_edist_tensor[i] = shortest_edist_tensor[i]-round(shortest_edist_tensor[i]/userInputs.domain_size[i])*userInputs.domain_size[i];
                            }
                        }
                        shortest_edist_tensor = userInputs.nucleation_parameters_list[userInputs.nucleation_parameters_list_index.at(thisNuclei->orderParameterIndex)].rotation_matrix * shortest_edist_tensor;
                        for (unsigned int i=0; i<dim; i++){
                            shortest_edist_tensor[i] /= userInputs.nucleation_parameters_list[userInputs.nucleation_parameters_list_index.at(thisNuclei->orderParameterIndex)].freeze_semiaxes[i];
                        }
                        weighted_dist = shortest_edist_tensor.norm_square();

						if (weighted_dist < 1.0 || thisNuclei->center.distance(q_point_list[q_point]) < diag_dist){
							if ((unsigned int)ti->level() < userInputs.max_refinement_level){
								mark_refine = true;
								break;
							}
						}
						if (mark_refine) break;
					}
					if (mark_refine) break;
				}
				if (mark_refine) di->set_refine_flag();
			}
			++di;
			++ti;
		}
		// The bulk of all of modifySolutionFields is spent in the following two function calls
		refineGrid();
		reinit();

		// If the mesh hasn't changed from the previous cycle, stop remeshing
		if (totalDOFs == numDoF_preremesh) break;
		numDoF_preremesh = totalDOFs;
	}
}

// Template instantiations
#include "../../include/matrixFreePDE_template_instantiations.h"
