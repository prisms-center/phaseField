// =================================================================================
// NUCLEATION FUNCTIONS
// =================================================================================

// =================================================================================
// Nucleation probability
// =================================================================================
template <int dim, int degree>
double customPDE<dim,degree>::nucProb(double cValue, double dV) const
{
	// Calculate the nucleation rate
	double J=k1*exp(-k2/(std::max(cValue,1.0e-6)) - tau/(this->currentTime));
	double retProb=1.0-exp(-J*timeStep*((double)skipNucleationSteps)*dV);
    return retProb;
}

// =================================================================================
// Get list of prospective new nuclei for the local processor
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::getLocalNucleiList(std::vector<nucleus<dim> > &newnuclei, std::vector<unsigned int> order_parameter_list, std::vector<unsigned int> other_var_list) const
{
	// Nickname for current time and time step
	double t=this->currentTime;
	unsigned int inc=this->currentIncrement;

    //QGauss<dim>  quadrature(degree+1);
    QGaussLobatto<dim>  quadrature(degree+1);
    FEValues<dim> fe_values (*(this->FESet[0]), quadrature, update_values|update_quadrature_points|update_JxW_values);
	const unsigned int   num_quad_points = quadrature.size();
	std::vector<std::vector<double> > op_values(order_parameter_list.size(),std::vector<double>(num_quad_points));
	std::vector<std::vector<double> > other_values(other_var_list.size(),std::vector<double>(num_quad_points));
	std::vector<dealii::Point<dim> > q_point_list(num_quad_points);

	std::vector<dealii::Point<dim> > q_point_list_overlap(num_quad_points);

	typename DoFHandler<dim>::active_cell_iterator   di = this->dofHandlersSet_nonconst[0]->begin_active();

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
            for (unsigned int var = 0; var < order_parameter_list.size(); var++){
            	fe_values.get_function_values(*(this->solutionSet[order_parameter_list[var]]), op_values[var]);
            }
            for (unsigned int var = 0; var < other_var_list.size(); var++){
            	fe_values.get_function_values(*(this->solutionSet[other_var_list[var]]), other_values[var]);
            }
            q_point_list = fe_values.get_quadrature_points();
            double ele_vol = 0.0;
            double ele_av_conc = 0.0;
        	// Loop over the quadrature points
            for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                ele_av_conc = ele_av_conc + other_values[0][q_point]*fe_values.JxW(q_point);
                ele_vol = ele_vol + fe_values.JxW(q_point);
            }
            ele_av_conc = ele_av_conc/ele_vol;
            
            //Compute random no. between 0 and 1 (new method)
            rand_val=distr(gen);
            //Nucleation probability
            double Prob=nucProb(ele_av_conc,ele_vol);

            //std::cout << "Nucleation probability: " << Prob << " Random value: " << rand_val << " Steady-state term: " << k1*exp(-k2/(std::max(ele_av_conc,1.0e-6))) << " Incubation term: " << k1*exp(- tau/(this->currentTime)) << std::endl;

        
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
                    bool periodic_j = (this->BC_list[1].var_BC_type[2*j]==PERIODIC);
                    bool insafetyzone_j = (periodic_j || ((nuc_ele_pos[j] > borderreg) && (nuc_ele_pos[j] < this->userInputs.domain_size[j]-borderreg)));
                    insafetyzone = insafetyzone && insafetyzone_j;
                }
                
                if (insafetyzone){
                
                	// Check to see if the order parameter anywhere within the element is above the threshold
                	bool anyqp_OK = false;
                	for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                		double sum_op = 0.0;
                		for (unsigned int num_op = 0; num_op < order_parameter_list.size(); num_op++){
                			sum_op += op_values[num_op][q_point];
                		}
                		if (sum_op < maxOrderParameterNucleation){
                			anyqp_OK =true;
                		}
                	}
                    
                	if (anyqp_OK){
                		// Pick the order parameter
                		std::random_device rd2;
                		std::mt19937 gen2(rd2());
                		std::uniform_int_distribution<unsigned int> int_distr(0,order_parameter_list.size()-1);
                		unsigned int op_for_nucleus = order_parameter_list[int_distr(gen2)];
                		std::cout << "Nucleation order parameter: " << op_for_nucleus << " " << rand_val << std::endl;
                		std::cout << "order parameter list: " << order_parameter_list[0] << " " << order_parameter_list[1] << " " << order_parameter_list[2] << std::endl;

                		//Add nucleus to prospective list
                		std::cout << "Prospective nucleation event. Nucleus no. " << nuclei.size()+1 << std::endl;
                		std::cout << "nucleus center " << nuc_ele_pos << std::endl;
                		nucleus<dim>* temp = new nucleus<dim>;
                		temp->index=nuclei.size();
                		temp->center=nuc_ele_pos;
                		temp->semiaxes.push_back(semiaxis_a);
                		temp->semiaxes.push_back(semiaxis_b);
                		if (dim == 3){
                			temp->semiaxes.push_back(semiaxis_c);
                		}
                		temp->seededTime=t;
                		temp->seedingTime = t_hold;
                		temp->seedingTimestep = inc;
                		temp->orderParameterIndex = op_for_nucleus;
                		newnuclei.push_back(*temp);
                	}
                }
            }
        }
    // Increment the cell iterators
    ++di;
	}
}

// =======================================================================================================
//Making sure all new nuclei from complete prospective list do not overlap with existing precipitates
// =======================================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::safetyCheckNewNuclei(std::vector<nucleus<dim> > newnuclei, std::vector<unsigned int> order_parameter_list, std::vector<unsigned int> &conflict_inds)
{
    //QGauss<dim>  quadrature(degree+1);
    QGaussLobatto<dim>  quadrature(degree+1);
    FEValues<dim> fe_values (*(this->FESet[0]), quadrature, update_values|update_quadrature_points|update_JxW_values);
    const unsigned int   num_quad_points = quadrature.size();
    std::vector<std::vector<double> > op_values(order_parameter_list.size(),std::vector<double>(num_quad_points));
    std::vector<dealii::Point<dim> > q_point_list(num_quad_points);

    //Nucleus cycle
    for (typename std::vector<nucleus<dim> >::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){
        bool isClose=false;

        //Element cycle
	    typename DoFHandler<dim>::active_cell_iterator   di = this->dofHandlersSet_nonconst[0]->begin_active();
        while (di != this->dofHandlersSet_nonconst[0]->end())
        {
            if (di->is_locally_owned()){
                fe_values.reinit(di);
                for (unsigned int var = 0; var < order_parameter_list.size(); var++){
                	fe_values.get_function_values(*(this->solutionSet[order_parameter_list[var]]), op_values[var]);
                }
                q_point_list = fe_values.get_quadrature_points();

                //Quadrature points cycle
                for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                    // Calculate the ellipsoidal distance to the center of the nucleus
                    double weighted_dist = 0.0;
                    for (unsigned int i=0; i<dim; i++){
                        double shortest_edist = thisNuclei->center(i) - q_point_list[q_point](i);
                        bool periodic_i = (this->BC_list[1].var_BC_type[2*i]==PERIODIC);
                        if (periodic_i){
                            double domsize =this->userInputs.domain_size[i];
                            shortest_edist = shortest_edist-round(shortest_edist/domsize)*domsize;
                        }
                        double temp = shortest_edist/(opfreeze_semiaxes[i]);
                        weighted_dist += temp*temp;
                    }
                    if (weighted_dist < 1.0){
                    	double sum_op = 0.0;
                    	for (unsigned int num_op = 0; num_op < order_parameter_list.size(); num_op++){
                    		sum_op += op_values[num_op][q_point];
                    	}
                        if (sum_op > 0.1){
                            isClose=true;
                            std::cout << "Attempted nucleation failed due to overlap w/ existing particle!!!!!!"  << std::endl;
                            conflict_inds.push_back(thisNuclei->index);
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
void customPDE<dim,degree>::refineMeshNearNuclei(std::vector<nucleus<dim> > newnuclei)
{
	//QGauss<dim>  quadrature(degree+1);
	QGaussLobatto<dim>  quadrature(degree+1);
	FEValues<dim> fe_values (*(this->FESet[0]), quadrature, update_values|update_quadrature_points|update_JxW_values);
	const unsigned int   num_quad_points = quadrature.size();
	std::vector<dealii::Point<dim> > q_point_list(num_quad_points);

	typename Triangulation<dim>::active_cell_iterator ti;
	typename DoFHandler<dim>::active_cell_iterator   di;

	unsigned int numDoF_preremesh = this->totalDOFs;

	for (unsigned int remesh_index=0; remesh_index < (this->userInputs.max_refinement_level-this->userInputs.min_refinement_level); remesh_index++){
		ti  = this->triangulation.begin_active();
		di = this->dofHandlersSet_nonconst[0]->begin_active();
		while (di != this->dofHandlersSet_nonconst[0]->end()){
			if (di->is_locally_owned()){

				bool mark_refine = false;

				fe_values.reinit (di);
				q_point_list = fe_values.get_quadrature_points();

				// Calculate the distance from the corner of the cell to the middle of the cell
				double diag_dist = 0.0;
				for (unsigned int i=0; i<dim; i++){
					diag_dist += (this->userInputs.domain_size[i]*this->userInputs.domain_size[i])/(this->userInputs.subdivisions[i]*this->userInputs.subdivisions[i]);
				}
				diag_dist = sqrt(diag_dist);
				diag_dist /= 2.0*pow(2.0,ti->level());

				for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
					for (typename std::vector<nucleus<dim> >::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){

						// Calculate the ellipsoidal distance to the center of the nucleus
						double weighted_dist = 0.0;
						for (unsigned int i=0; i<dim; i++){
							double shortest_edist = thisNuclei->center(i) - q_point_list[q_point](i);
							bool periodic_i = (this->BC_list[1].var_BC_type[2*i]==PERIODIC);
							if (periodic_i){
								double domsize =this->userInputs.domain_size[i];
								shortest_edist = shortest_edist-round(shortest_edist/domsize)*domsize;
							}
							double temp = shortest_edist/(opfreeze_semiaxes[i]);
							weighted_dist += temp*temp;
						}

						if (weighted_dist < 1.0 || thisNuclei->center.distance(q_point_list[q_point]) < diag_dist){
							if ((unsigned int)ti->level() < this->userInputs.max_refinement_level){
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
		this->refineGrid();
		this->reinit();

		// If the mesh hasn't changed from the previous cycle, stop remeshing
		if (this->totalDOFs == numDoF_preremesh) break;
		numDoF_preremesh = this->totalDOFs;
	}
}

// =================================================================================
// Global nucleation procedure
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::getNucleiList()
{
    if ( this->currentIncrement % skipNucleationSteps == 0){

		// Declare a vector of all the NEW nuclei seeded in this time step
		std::vector<nucleus<dim> > newnuclei;

		// Get list of prospective new nuclei for the local processor
		this->pcout << "Nucleation attempt for increment " << this->currentIncrement << std::endl;
		std::vector<unsigned int> order_parameter_list;
		std::vector<unsigned int> other_var_list;
		other_var_list.push_back(0);
		order_parameter_list.push_back(1);
		order_parameter_list.push_back(2);
		order_parameter_list.push_back(3);
		getLocalNucleiList(newnuclei,order_parameter_list,other_var_list);

		// Generate global list of new nuclei and resolve conflicts between new nuclei
		parallelNucleationList<dim> new_nuclei_parallel(newnuclei);
		newnuclei = new_nuclei_parallel.buildGlobalNucleiList(minDistBetweenNuclei, nuclei.size());

		// Final check to resolve overlap conflicts with existing precipitates
		std::vector<unsigned int> conflict_inds;
		safetyCheckNewNuclei(newnuclei, order_parameter_list, conflict_inds);

		newnuclei = new_nuclei_parallel.removeSubsetOfNuclei(conflict_inds, nuclei.size());

        // Add the new nuclei to the list of nuclei
        nuclei.insert(nuclei.end(),newnuclei.begin(),newnuclei.end());

        // Refine mesh near the new nuclei
        if (newnuclei.size() > 0 && this->userInputs.h_adaptivity == true){
        	refineMeshNearNuclei(newnuclei);
        }
    }
    else if (this->currentIncrement == 1){
    	// Declare a vector of all the NEW nuclei seeded in this time step
    	std::vector<nucleus<dim> > newnuclei;

    	// Get list of prospective new nuclei for the local processor
    	this->pcout << "Nucleation attempt for increment " << this->currentIncrement << std::endl;
    	std::vector<unsigned int> order_parameter_list;
    	std::vector<unsigned int> other_var_list;
    	other_var_list.push_back(0);
    	order_parameter_list.push_back(1);
    	order_parameter_list.push_back(2);
    	order_parameter_list.push_back(3);

    	if (this->currentIncrement == 1){
    		this->pcout << "newnuclei size " << newnuclei.size() << std::endl;
    		while (newnuclei.size() == 0){
    			this->currentTime+=this->userInputs.dtValue*(double)skipNucleationSteps;
    			this->currentIncrement+=skipNucleationSteps;

    			while (this->outputTimeStepList.size() > 0 && this->outputTimeStepList[0] < this->currentIncrement){
    				this->outputTimeStepList.erase(this->outputTimeStepList.begin());
    			}

    			getLocalNucleiList(newnuclei,order_parameter_list,other_var_list);
    			this->pcout << "nucleation attempt! " << this->currentTime << " " << this->currentIncrement << std::endl;


    			// Generate global list of new nuclei and resolve conflicts between new nuclei
    			parallelNucleationList<dim> new_nuclei_parallel(newnuclei);
    			newnuclei = new_nuclei_parallel.buildGlobalNucleiList(minDistBetweenNuclei, nuclei.size());

    			// Final check to resolve overlap conflicts with existing precipitates
    			std::vector<unsigned int> conflict_inds;
    			safetyCheckNewNuclei(newnuclei, order_parameter_list, conflict_inds);

    			newnuclei = new_nuclei_parallel.removeSubsetOfNuclei(conflict_inds, nuclei.size());

    		}
    	}

    	// Add the new nuclei to the list of nuclei
    	nuclei.insert(nuclei.end(),newnuclei.begin(),newnuclei.end());

    	// Refine mesh near the new nuclei
    	if (newnuclei.size() > 0 && this->userInputs.h_adaptivity == true){
    		refineMeshNearNuclei(newnuclei);
    	}
    }

}
