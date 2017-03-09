// =================================================================================
// NUCLEATION FUNCTIONS
// =================================================================================


// =================================================================================
// Nucleation probability
// =================================================================================
template <int dim, int degree>
double customPDE<dim,degree>::nucProb(double cValue, double dV) const
{
	//Nucleation rate
	double J=k1*exp(-k2/(std::max(cValue-calmin,1.0e-6)));
	//We need element volume (or area in 2D)
	double retProb=1.0-exp(-J*timeStep*((double)skipNucleationSteps)*dV);
    return retProb;
}


// =================================================================================
// Sends the list of new nuclei to the next processor
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::sendUpdate (std::vector<nucleus> &newnuclei, int procno) const
{
    int currnonucs=newnuclei.size();
    //MPI SECTION TO SEND INFORMATION TO THE PROCESSOR procno
    //Sending local no. of nuclei
    MPI_Send(&currnonucs, 1, MPI_INT, procno, 0, MPI_COMM_WORLD);
    if (currnonucs > 0){
        //Creating vectors of each quantity in nuclei. Each numbered acording to the tags used for MPI_Send/MPI_Recv
        //1 - index
        std::vector<unsigned int> s_index;
        //2 - "x" componenet of center
        std::vector<double> s_center_x;
        //3 - "y" componenet of center
        std::vector<double> s_center_y;
        //4 - "z" componenet of center
        std::vector<double> s_center_z;
        //5 - radius
        std::vector<double> s_radius;
        //6 - seededTime
        std::vector<double> s_seededTime;
        //7 - seedingTime
        std::vector<double> s_seedingTime;
        //8 - seedingTimestep
        std::vector<unsigned int> s_seedingTimestep;
        
        //Loop to store info of all nuclei into vectors
        for (std::vector<nucleus>::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){
            s_index.push_back(thisNuclei->index);
            dealii::Point<problemDIM> s_center=thisNuclei->center;
            s_center_x.push_back(s_center[0]);
            s_center_y.push_back(s_center[1]);
            s_center_z.push_back(s_center[2]);
            s_radius.push_back(thisNuclei->radius);
            s_seededTime.push_back(thisNuclei->seededTime);
            s_seedingTime.push_back(thisNuclei->seedingTime);
            s_seedingTimestep.push_back(thisNuclei->seedingTimestep);
        }
        //Send vectors to next processor
        MPI_Send(&s_index[0], currnonucs, MPI_UNSIGNED, procno, 1, MPI_COMM_WORLD);
        MPI_Send(&s_center_x[0], currnonucs, MPI_DOUBLE, procno, 2, MPI_COMM_WORLD);
        MPI_Send(&s_center_y[0], currnonucs, MPI_DOUBLE, procno, 3, MPI_COMM_WORLD);
        MPI_Send(&s_center_z[0], currnonucs, MPI_DOUBLE, procno, 4, MPI_COMM_WORLD);
        MPI_Send(&s_radius[0], currnonucs, MPI_DOUBLE, procno, 5, MPI_COMM_WORLD);
        MPI_Send(&s_seededTime[0], currnonucs, MPI_DOUBLE, procno, 6, MPI_COMM_WORLD);
        MPI_Send(&s_seedingTime[0], currnonucs, MPI_DOUBLE, procno, 7, MPI_COMM_WORLD);
        MPI_Send(&s_seedingTimestep[0], currnonucs, MPI_UNSIGNED, procno, 8, MPI_COMM_WORLD);
    }
    //END OF MPI SECTION
}

// =================================================================================
// Recieves the list of new nuclei to the next processor
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::receiveUpdate (std::vector<nucleus> &newnuclei, int procno) const
{
    //MPI PROCEDURE TO RECIEVE INFORMATION FROM ANOTHER PROCESSOR AND UPDATE LOCAL NUCLEI INFORMATION
    int recvnonucs = 0;
    int currnonucs = newnuclei.size();
    MPI_Recv(&recvnonucs, 1, MPI_INT, procno, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (recvnonucs > 0){

        //Creating vectors of each quantity in nuclei. Each numbered acording to the tags used for MPI_Send/MPI_Recv
        //1 - index
        std::vector<unsigned int> r_index(recvnonucs,0);
        //2 - "x" componenet of center
        std::vector<double> r_center_x(recvnonucs,0.0);
        //3 - "y" componenet of center
        std::vector<double> r_center_y(recvnonucs,0.0);
        //4 - "z" componenet of center
        std::vector<double> r_center_z(recvnonucs,0.0);
        //5 - radius
        std::vector<double> r_radius(recvnonucs,0.0);
        //6 - seededTime
        std::vector<double> r_seededTime(recvnonucs,0.0);
        //7 - seedingTime
        std::vector<double> r_seedingTime(recvnonucs,0.0);
        //8 - seedingTimestep
        std::vector<unsigned int> r_seedingTimestep(recvnonucs,0);

        //Recieve vectors from processor procno
        MPI_Recv(&r_index[0], recvnonucs, MPI_UNSIGNED, procno, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_center_x[0], recvnonucs, MPI_DOUBLE, procno, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_center_y[0], recvnonucs, MPI_DOUBLE, procno, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_center_z[0], recvnonucs, MPI_DOUBLE, procno, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_radius[0], recvnonucs, MPI_DOUBLE, procno, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_seededTime[0], recvnonucs, MPI_DOUBLE, procno, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_seedingTime[0], recvnonucs, MPI_DOUBLE, procno, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_seedingTimestep[0], recvnonucs, MPI_UNSIGNED, procno, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //Loop to store info in vectors onto the nuclei structure
        for (int jnuc=0; jnuc<=recvnonucs-1; jnuc++){
            nucleus* temp = new nucleus;
            temp->index=r_index[jnuc];
            dealii::Point<problemDIM> r_center;
            r_center[0]=r_center_x[jnuc];
            r_center[1]=r_center_y[jnuc];
            r_center[2]=r_center_z[jnuc];
            temp->center=r_center;
            temp->radius=r_radius[jnuc];
            temp->seededTime=r_seededTime[jnuc];
            temp->seedingTime = r_seedingTime[jnuc];
            temp->seedingTimestep = r_seedingTimestep[jnuc];
            newnuclei.push_back(*temp);
        }

    }
}

// =================================================================================
// Broadcast the final list of new nuclei from the last processor to the rest
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::broadcastUpdate (std::vector<nucleus> &newnuclei, int broadcastProc, int thisProc) const
{
    //MPI PROCEDURE TO SEND THE LIST OF NEW NUCLEI FROM ONE PROCESSOR TO ALL THE OTHERS
    int currnonucs = newnuclei.size();
    MPI_Bcast(&currnonucs, 1, MPI_INT, broadcastProc, MPI_COMM_WORLD);
    if (currnonucs > 0){
        
        //Creating vectors of each quantity in nuclei. Each numbered acording to the tags used for MPI_Send/MPI_Recv
    	unsigned int initial_vec_size;
    	if (thisProc == broadcastProc){
    		initial_vec_size = 0;
    	}
    	else{
    		initial_vec_size = currnonucs;
    	}

        //1 - index
        std::vector<unsigned int> r_index(initial_vec_size,0);
        //2 - "x" componenet of center
        std::vector<double> r_center_x(initial_vec_size,0.0);
        //3 - "y" componenet of center
        std::vector<double> r_center_y(initial_vec_size,0.0);
        //4 - "z" componenet of center
        std::vector<double> r_center_z(initial_vec_size,0.0);
        //5 - radius
        std::vector<double> r_radius(initial_vec_size,0.0);
        //6 - seededTime
        std::vector<double> r_seededTime(initial_vec_size,0.0);
        //7 - seedingTime
        std::vector<double> r_seedingTime(initial_vec_size,0.0);
        //8 - seedingTimestep
        std::vector<unsigned int> r_seedingTimestep(initial_vec_size,0);
        
        if (thisProc == broadcastProc){
        	for (std::vector<nucleus>::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){
        		r_index.push_back(thisNuclei->index);
        		dealii::Point<problemDIM> s_center=thisNuclei->center;
        		r_center_x.push_back(s_center[0]);
        		r_center_y.push_back(s_center[1]);
                r_center_z.push_back(s_center[2]);
        		r_radius.push_back(thisNuclei->radius);
        		r_seededTime.push_back(thisNuclei->seededTime);
        		r_seedingTime.push_back(thisNuclei->seedingTime);
        		r_seedingTimestep.push_back(thisNuclei->seedingTimestep);
        	}
        }

        //Recieve vectors from processor procno
        MPI_Bcast(&r_index[0], currnonucs, MPI_UNSIGNED, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_center_x[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_center_y[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_center_z[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_radius[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_seededTime[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_seedingTime[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_seedingTimestep[0], currnonucs, MPI_UNSIGNED, broadcastProc, MPI_COMM_WORLD);
        
        newnuclei.clear();

        //Loop to store info in vectors onto the nuclei structure
        for (int jnuc=0; jnuc<=currnonucs-1; jnuc++){
            nucleus* temp = new nucleus;
            temp->index=r_index[jnuc];
            dealii::Point<problemDIM> r_center;
            r_center[0]=r_center_x[jnuc];
            r_center[1]=r_center_y[jnuc];
            r_center[2]=r_center_z[jnuc];
            temp->center=r_center;
            temp->radius=r_radius[jnuc];
            temp->seededTime=r_seededTime[jnuc];
            temp->seedingTime = r_seedingTime[jnuc];
            temp->seedingTimestep = r_seedingTimestep[jnuc];
            newnuclei.push_back(*temp);
        }
        
    }
}

// =================================================================================
// Determine if any new nuclei are in conflict and resolve those conflicts
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::resolveNucleationConflicts (std::vector<nucleus> &newnuclei) const{

	std::vector<nucleus> newnuclei_cleaned;
	unsigned int old_num_nuclei = nuclei.size();

	for (unsigned int nuc_index=0; nuc_index<newnuclei.size(); nuc_index++){
		bool isClose=false;

		for (unsigned int prev_nuc_index=0; prev_nuc_index<nuc_index; prev_nuc_index++){

			// We may want to break this section into a separate function to allow different choices for when
			// nucleation should be prevented
			if (newnuclei[nuc_index].center.distance(newnuclei[prev_nuc_index].center) < minDistBetwenNuclei){
				isClose = true;
				std::cout << "Conflict between nuclei! Distance is: " << newnuclei[nuc_index].center.distance(newnuclei[prev_nuc_index].center)
						<< " Conflict removed."<< std::endl;
				break;
			}
		}

		if (!isClose){
			newnuclei[nuc_index].index = old_num_nuclei + newnuclei_cleaned.size();
			newnuclei_cleaned.push_back(newnuclei[nuc_index]);
		}
	}

	newnuclei = newnuclei_cleaned;
}

// =================================================================================
// Get list of prospective new nuclei for the local processor
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::getLocalNucleiList(std::vector<nucleus> &newnuclei) const
{
	// Nickname for current time and time step
	double t=this->currentTime;
	unsigned int inc=this->currentIncrement;


	QGauss<dim>  quadrature(degree+1);
	FEValues<dim> fe_values (*(this->FESet[0]), quadrature, update_values|update_quadrature_points|update_JxW_values);
	const unsigned int   num_quad_points = quadrature.size();
	std::vector<double> var_value(num_quad_points);
	std::vector<double> var_value2(num_quad_points);
	std::vector<dealii::Point<dim> > q_point_list(num_quad_points);

	std::vector<dealii::Point<dim> > q_point_list_overlap(num_quad_points);
	std::vector<double> var_value2_overlap(num_quad_points);

	typename DoFHandler<dim>::active_cell_iterator   di = this->dofHandlersSet_nonconst[0]->begin_active();

	// What used to be in nuc_attempt
	double rand_val;
	//Better random no. generator
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> distr(0.0,1.0);

	while (di != this->dofHandlersSet_nonconst[0]->end())
	{
		if (di->is_locally_owned()){

			fe_values.reinit (di);
			fe_values.get_function_values(*(this->solutionSet[0]), var_value);
			fe_values.get_function_values(*(this->solutionSet[1]), var_value2);
			q_point_list = fe_values.get_quadrature_points();

			// Loop over the quadrature points
			for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
				bool insafetyzone = (q_point_list[q_point][0] > borderreg) && (q_point_list[q_point][0] < spanX-borderreg) && (q_point_list[q_point][1] > borderreg) && (q_point_list[q_point][1] < spanY-borderreg) && (q_point_list[q_point][2] > borderreg) && (q_point_list[q_point][2] < spanY-borderreg);
				bool periodic=false;
				if (insafetyzone || periodic){
					if (var_value2[q_point] < maxOrderParameterNucleation){

						//Compute random no. between 0 and 1 (old method)
						rand_val=distr(gen);
						double Prob=nucProb(var_value[q_point],fe_values.JxW(q_point));
						if (rand_val <= Prob){
							//std::cout << "random value " << rand_val << ", probability " << Prob << std::endl;

							//loop over all existing nuclei to check if they are in the vicinity
							bool isClose=false;

							// ---------------------------------------------------------------------------------
							// Check to see if the prospective nucleus would overlap with any other particles
							// Not necessarily needed, but good to keep around
	//			    			typename DoFHandler<dim>::active_cell_iterator   di_overlap = this->dofHandlersSet_nonconst[0]->begin_active();
	//			    			while (di_overlap != this->dofHandlersSet_nonconst[0]->end())
	//			    			{
	//			    				if (di_overlap->is_locally_owned()){
	//			    					fe_values.reinit (di_overlap);
	//			    					fe_values.get_function_values(*(this->solutionSet[1]), var_value2_overlap);
	//			    					q_point_list_overlap = fe_values.get_quadrature_points();
	//			    					for (unsigned int q_point_overlap=0; q_point_overlap<num_quad_points; ++q_point_overlap){
	//			    						if (q_point_list_overlap[q_point_overlap].distance(q_point_list[q_point])<opfreeze_radius){
	//			    							if (var_value2_overlap[q_point_overlap] > 0.1){
	//			    								isClose=true;
	//			    								std::cout << "Attempted nucleation failed due to overlap w/ existing particle!!!!!!" << q_point_list_overlap[q_point_overlap](0) << " " << q_point_list_overlap[q_point_overlap](1) << " " << q_point_list[q_point](0) << " " << q_point_list[q_point](1) << std::endl;
	//			    								break;
	//			    							}
	//			    						}
	//			    					}
	//			    					if (isClose) break;
	//			    				}
	//			    				++di_overlap;
	//			    			}
	//			    			fe_values.reinit (di);
							// ---------------------------------------------------------------------------------

							if (!isClose){
								std::cout << "Nucleation event. Nucleus no. " << nuclei.size()+1 << std::endl;
								std::cout << "nucleus center " << q_point_list[q_point] << std::endl;
								nucleus* temp = new nucleus;
								temp->index=nuclei.size();
								temp->center=q_point_list[q_point];
								temp->radius=n_radius;
								temp->seededTime=t;
								temp->seedingTime = t_hold;
								temp->seedingTimestep = inc;
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

// =================================================================================
// Generate global list of new nuclei and resolve conflicts between new nuclei
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::buildGlobalNucleiList(std::vector<nucleus> &newnuclei) const
{
	//MPI INITIALIZATON
	int numProcs=Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
	int thisProc=Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
	if (numProcs > 1) {
		// Cycle through each processor, sending and receiving, to append the list of new nuclei
		for (int proc_index=0; proc_index < numProcs-1; proc_index++){
			if (thisProc == proc_index){
				sendUpdate(newnuclei, thisProc+1);
			}
			else if (thisProc == proc_index+1){
				receiveUpdate(newnuclei, thisProc-1);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		// The final processor now has all of the new nucleation attempts
		// Check for conflicts on the final processor before broadcasting the list
		if (thisProc == numProcs-1){
			resolveNucleationConflicts(newnuclei);
		}

		// The final processor now has the final list of the new nuclei, broadcast it to all the other processors
		broadcastUpdate(newnuclei, numProcs-1, thisProc);
	}
	else {
		// Check for conflicts between nucleation attempts this time step
		resolveNucleationConflicts(newnuclei);
	}
}

// =================================================================================
// Refine mesh near the new nuclei
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::refineMeshNearNuclei(std::vector<nucleus> newnuclei)
{
	QGauss<dim>  quadrature(degree+1);
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

				for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
					for (std::vector<nucleus>::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){
						if (thisNuclei->center.distance(q_point_list[q_point])<opfreeze_radius){
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
    if ( this->currentIncrement % skipNucleationSteps == 0 ){

		// Declare alize vector of all the NEW nuclei seeded in this time step
		std::vector<nucleus> newnuclei;

		// Get list of prospective new nuclei for the local processor
		std::cout << "Nucleation attempt for increment " << this->currentIncrement << std::endl;
		getLocalNucleiList(newnuclei);

		// Generate global list of new nuclei and resolve conflicts between new nuclei
		buildGlobalNucleiList(newnuclei);

        // Add the new nuclei to the list of nuclei
        nuclei.insert(nuclei.end(),newnuclei.begin(),newnuclei.end());

        // Refine mesh near the new nuclei
        if (newnuclei.size() > 0 && this->userInputs.h_adaptivity == true){
        	refineMeshNearNuclei(newnuclei);
        }
    }
}
