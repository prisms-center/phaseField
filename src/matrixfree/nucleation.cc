// Methods in MatrixFreePDE to update the list of nuclei
#include "../../include/matrixFreePDE.h"
#include "../../include/parallelNucleationList.h"
#include "../../include/varBCs.h"
#include <random>
#include <time.h>

// =======================================================================================================
// Function called in solve to update the global list of nuclei
// =======================================================================================================
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::updateNucleiList() {

  if (userInputs.nucleation_occurs){
      
      
    //
    if (currentIncrement % userInputs.steps_between_nucleation_attempts == 0 || currentIncrement == 1){
      
      if (userInputs.dtValue*(double)currentIncrement >= userInputs.nucleation_start_time && userInputs.dtValue*(double)currentIncrement <= userInputs.nucleation_end_time){
        computing_timer.enter_subsection("matrixFreePDE: nucleation");
        
        // Apply constraints
        for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
          constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
          constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
          solutionSet[fieldIndex]->update_ghost_values();
        }
        
        std::vector<nucleus<dim> > new_nuclei;
        
        // Disabling skipping steps
        /*
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
         */
        //NEW line below instead of if statement (uncomment later)
        new_nuclei = getNewNuclei();
        nuclei.insert(nuclei.end(),new_nuclei.begin(),new_nuclei.end());
       // std::cout << "Total Nucleus no. " << (new_nuclei.size()) << std::endl;
        if (new_nuclei.size() > 0 && userInputs.h_adaptivity == true){
          refineMeshNearNuclei(new_nuclei);
        }
        
        computing_timer.leave_subsection("matrixFreePDE: nucleation");
      }
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
    newnuclei = new_nuclei_parallel.buildGlobalNucleiList(userInputs.min_distance_between_nuclei, userInputs.min_distance_between_OP, nuclei.size());

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
            // NOTE: This might not be the best way to do this. This is missing the loop of the DoFs from Step-3

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
            
                // NEW SECTION: Loop through each existing nuclei to verify if any of them correspond to the current order parameter
                unsigned int numbercurrnuclei = 0;
                for (unsigned int nucind = 0; nucind < nuclei.size(); nucind++){
                    if (nuclei[nucind].orderParameterIndex==variable_index){
                        numbercurrnuclei++;
                        //pcout << "There are already " << numbercurrnuclei << " for order parameter " << variable_index << std::endl;
                    }
                }
                // If only one nucleus per order parameter is allowed do not enter the next section
                if ((numbercurrnuclei == 0) || (userInputs.multiple_nuclei_per_order_parameter == true)){
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
                        for (unsigned int j=0; j<dim; j++){
                            randvec[j]=distr(gen);
                            nuc_ele_pos[j]=ele_origin[j] + (ele_max[j]-ele_origin[j])*randvec[j];
                        }

                        //Make sure point is in safety zone
                        bool insafetyzone = true;
                        for (unsigned int j=0; j < dim; j++){
                            bool periodic_j = (userInputs.BC_list[1].var_BC_type[2*j]==PERIODIC);
                            bool insafetyzone_j = (periodic_j || ((nuc_ele_pos[j] > userInputs.get_no_nucleation_border_thickness(variable_index)) && (nuc_ele_pos[j] < userInputs.domain_size[j]-userInputs.get_no_nucleation_border_thickness(variable_index))));
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
                                temp->semiaxes = userInputs.get_nucleus_semiaxes(variable_index);
                                temp->seededTime=t;
                                temp->seedingTime = userInputs.get_nucleus_hold_time(variable_index);
                                temp->seedingTimestep = inc;
                                temp->orderParameterIndex = variable_index;
                                newnuclei.push_back(*temp);
                            }
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

    //NEW section check if order parameters from prospective nuclei overlap with those from existing nuclei
    if (userInputs.multiple_nuclei_per_order_parameter == false){
      bool opOverlap;
        for (typename std::vector<nucleus<dim> >::iterator thisNucleus1=newnuclei.begin(); thisNucleus1!=newnuclei.end(); ++thisNucleus1){
            opOverlap = false;
            for (typename std::vector<nucleus<dim> >::iterator thisNucleus2=thisNucleus1+1; thisNucleus2!=newnuclei.end(); ++thisNucleus2){
                if(thisNucleus1->orderParameterIndex == thisNucleus2->orderParameterIndex){
                    opOverlap=true;
                    std::cout << "Attempted nucleation failed due to overlap w/ existing order parameter!"  << std::endl;
                    conflict_ids.push_back(thisNucleus1->index);
                }
            if (opOverlap) break;
            }
        }
    }
    // End of overlap of order parameters check
  
    //Nucleus cycle to check for overlapping nuclei
    for (typename std::vector<nucleus<dim> >::iterator thisNucleus=newnuclei.begin(); thisNucleus!=newnuclei.end(); ++thisNucleus){
        bool isClose=false;

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
                    double weighted_dist = weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_freeze_semiaxes(thisNucleus->orderParameterIndex), q_point_list[q_point], thisNucleus->orderParameterIndex);

                    if (weighted_dist < 1.0){
                    	double sum_op = 0.0;
                    	for (unsigned int num_op = 0; num_op < userInputs.nucleating_variable_indices.size(); num_op++){
                    		sum_op += op_values[num_op][q_point];
                    	}
                        if (sum_op > 0.1){
                            isClose=true;
                            std::cout << "Attempted nucleation failed due to overlap w/ existing particle!!!!!!"  << std::endl;
                            conflict_ids.push_back(thisNucleus->index);
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
					for (typename std::vector<nucleus<dim> >::iterator thisNucleus=newnuclei.begin(); thisNucleus!=newnuclei.end(); ++thisNucleus){

                        // Calculate the ellipsoidal distance to the center of the nucleus
                        double weighted_dist = weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_freeze_semiaxes(thisNucleus->orderParameterIndex), q_point_list[q_point], thisNucleus->orderParameterIndex);

						if (weighted_dist < 1.0 || thisNucleus->center.distance(q_point_list[q_point]) < diag_dist){
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

// First of two versions of this function to calculated the weighted distance from the center of a nucleus_semiaxes
// This version is for when the points are given as doubles
template <int dim, int degree>
double MatrixFreePDE<dim,degree>::weightedDistanceFromNucleusCenter(
    const dealii::Point<dim,double> center,
    const std::vector<double> semiaxes,
    const dealii::Point<dim,double> q_point_loc,
    const unsigned int var_index) const
{
    double weighted_dist = 0.0;
    dealii::Tensor<1,dim,double> shortest_edist_tensor = center - q_point_loc;
    for (unsigned int i=0; i<dim; i++){
        if (userInputs.BC_list[var_index].var_BC_type[2*i]==PERIODIC){
            shortest_edist_tensor[i] = shortest_edist_tensor[i]-round(shortest_edist_tensor[i]/userInputs.domain_size[i])*userInputs.domain_size[i];
        }
    }
    shortest_edist_tensor = userInputs.get_nucleus_rotation_matrix(var_index) * shortest_edist_tensor;
    for (unsigned int i=0; i<dim; i++){
        shortest_edist_tensor[i] /= semiaxes[i];
    }
    weighted_dist = shortest_edist_tensor.norm();
    return weighted_dist;
}

// Second of two versions of this function to calculated the weighted distance from the center of a nucleus_semiaxes
// This version is for when the points are given as vectorized arrays
template <int dim, int degree>
dealii::VectorizedArray<double> MatrixFreePDE<dim,degree>::weightedDistanceFromNucleusCenter(
    const dealii::Point<dim,double> center,
    const std::vector<double> semiaxes,
    const dealii::Point<dim,dealii::VectorizedArray<double> > q_point_loc,
    const unsigned int var_index) const
{
    dealii::VectorizedArray<double> weighted_dist = constV(0.0);
    dealii::Tensor<1,dim,dealii::VectorizedArray<double> > shortest_edist_tensor;
    for (unsigned int j=0; j<dim; j++){
        shortest_edist_tensor[j] = center(j) - q_point_loc(j); // Can I do this outside the loop?

        if (userInputs.BC_list[var_index].var_BC_type[2*j]==PERIODIC){
            for (unsigned k=0; k<q_point_loc(0).size();k++){
                shortest_edist_tensor[j][k] = shortest_edist_tensor[j][k]-round(shortest_edist_tensor[j][k]/userInputs.domain_size[j])*userInputs.domain_size[j];
            }
        }
    }
    shortest_edist_tensor = userInputs.get_nucleus_rotation_matrix(var_index) * shortest_edist_tensor;
    for (unsigned int j=0; j<dim; j++){
        shortest_edist_tensor[j] /= constV(semiaxes[j]);
    }
    weighted_dist = shortest_edist_tensor.norm_square();
    for (unsigned k=0; k<q_point_loc(0).size();k++){
        weighted_dist[k] = sqrt(weighted_dist[k]);
    }
    return weighted_dist;
}

// Template instantiations
#include "../../include/matrixFreePDE_template_instantiations.h"
