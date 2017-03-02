// =================================================================================
//Function prototypes
// =================================================================================
double nucProb(double cValue);
void nucAttempt(std::vector<nucleus> &newnuclei, std::map<dealii::types::global_dof_index, dealii::Point<problemDIM> > support_points, vectorType* c, vectorType* n, double t, unsigned int inc);
void receiveUpdate (std::vector<nucleus> &newnuclei, int procno);
void sendUpdate (std::vector<nucleus> &newnuclei, int procno);
void broadcastUpdate (std::vector<nucleus> &newnuclei, int broadcastProc, int thisProc);
void resolveNucleationConflicts(std::vector<nucleus> &newnuclei);

// =================================================================================
// NUCLEATION FUNCTIONS
// =================================================================================

// =================================================================================
// Nucleation probability
// =================================================================================
double nucProb(double cValue)
{
    //minimum grid spacing
    double dx=spanX/(std::pow(2.0,refineFactor)*(double)finiteElementDegree);
	//Nucleation rate
	double J=k1*exp(-k2/(cValue-calmin))*10.0;
	//We need element volume (or area in 2D)
	double retProb=1.0-exp(-J*timeStep*((double)skipNucleationSteps)*dx*dx);
    return retProb;
}

// =================================================================================
// Nucleation Attempt
// =================================================================================
void nucAttempt(std::vector<nucleus> &newnuclei, std::map<dealii::types::global_dof_index, dealii::Point<problemDIM> > support_points, vectorType* c, vectorType* n, double t, unsigned int inc)
{
    double rand_val;
    std::cout << "nucleation attempt" << std::endl;
    //Better random no. generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distr(0.0,1.0);
    //add nuclei based on concentration field values
    //loop over all points in the domain
    for (typename std::map<dealii::types::global_dof_index, dealii::Point<problemDIM> >::iterator it=support_points.begin(); it!=support_points.end(); ++it){
        unsigned int dof=it->first;
        //set only local owned values of the parallel vector
        if (n->locally_owned_elements().is_element(dof)){
            dealii::Point<problemDIM> nodePoint=it->second;
            //Safety zone to avoid nucleation near the edges in no-flux BCs
            bool insafetyzone = (nodePoint[0] > borderreg) && (nodePoint[0] < spanX-borderreg) && (nodePoint[1] > borderreg) && (nodePoint[1] < spanY-borderreg);
            //bool periodic = ((BC_list[0].var_BC_type == "PERIODIC") && (BC_list[1].var_BC_type == "PERIODIC"));
            bool periodic=false;
            if (insafetyzone || periodic){
                double nValue=(*n)(dof);
                double cValue=(*c)(dof);
                //Compute random no. between 0 and 1 (old method)
                rand_val=distr(gen);
                //Nucleation probability for element
                double Prob=nucProb(cValue);
                if (rand_val <= Prob){
                    //std::cout << "random value " << rand_val << ", probability " << Prob << std::endl;
                    //loop over all existing nuclei to check if they are in the vicinity
                    bool isClose=false;
                    for (std::vector<nucleus>::iterator thisNuclei=nuclei.begin(); thisNuclei!=nuclei.end(); ++thisNuclei){
                        if (thisNuclei->center.distance(nodePoint)<minDistBetwenNuclei){
                            isClose=true;
                        }
                    }
                    if (!isClose){
                        std::cout << "Nucleation event. Nucleus no. " << nuclei.size()+1 << std::endl;
                        std::cout << "nucleus center " << nodePoint << std::endl;
                        nucleus* temp = new nucleus;
                        temp->index=nuclei.size();
                        temp->center=nodePoint;
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

void sendUpdate (std::vector<nucleus> &newnuclei, int procno)
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
        //4 - radius
        std::vector<double> s_radius;
        //5 - seededTime
        std::vector<double> s_seededTime;
        //6 - seedingTime
        std::vector<double> s_seedingTime;
        //7 - seedingTimestep
        std::vector<unsigned int> s_seedingTimestep;
        
        //Loop to store info of all nuclei into vectors
        for (std::vector<nucleus>::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){
            s_index.push_back(thisNuclei->index);
            dealii::Point<problemDIM> s_center=thisNuclei->center;
            s_center_x.push_back(s_center[0]);
            s_center_y.push_back(s_center[1]);
            s_radius.push_back(thisNuclei->radius);
            s_seededTime.push_back(thisNuclei->seededTime);
            s_seedingTime.push_back(thisNuclei->seedingTime);
            s_seedingTimestep.push_back(thisNuclei->seedingTimestep);
        }
        //Send vectors to next processor
        MPI_Send(&s_index[0], currnonucs, MPI_UNSIGNED, procno, 1, MPI_COMM_WORLD);
        MPI_Send(&s_center_x[0], currnonucs, MPI_DOUBLE, procno, 2, MPI_COMM_WORLD);
        MPI_Send(&s_center_y[0], currnonucs, MPI_DOUBLE, procno, 3, MPI_COMM_WORLD);
        MPI_Send(&s_radius[0], currnonucs, MPI_DOUBLE, procno, 4, MPI_COMM_WORLD);
        MPI_Send(&s_seededTime[0], currnonucs, MPI_DOUBLE, procno, 5, MPI_COMM_WORLD);
        MPI_Send(&s_seedingTime[0], currnonucs, MPI_DOUBLE, procno, 6, MPI_COMM_WORLD);
        MPI_Send(&s_seedingTimestep[0], currnonucs, MPI_UNSIGNED, procno, 7, MPI_COMM_WORLD);
    }
    //END OF MPI SECTION
}

void receiveUpdate (std::vector<nucleus> &newnuclei, int procno)
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
        //4 - radius
        std::vector<double> r_radius(recvnonucs,0.0);
        //5 - seededTime
        std::vector<double> r_seededTime(recvnonucs,0.0);
        //6 - seedingTime
        std::vector<double> r_seedingTime(recvnonucs,0.0);
        //7 - seedingTimestep
        std::vector<unsigned int> r_seedingTimestep(recvnonucs,0);

        //Recieve vectors from processor procno
        MPI_Recv(&r_index[0], recvnonucs, MPI_UNSIGNED, procno, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_center_x[0], recvnonucs, MPI_DOUBLE, procno, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_center_y[0], recvnonucs, MPI_DOUBLE, procno, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_radius[0], recvnonucs, MPI_DOUBLE, procno, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_seededTime[0], recvnonucs, MPI_DOUBLE, procno, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_seedingTime[0], recvnonucs, MPI_DOUBLE, procno, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_seedingTimestep[0], recvnonucs, MPI_UNSIGNED, procno, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //Loop to store info in vectors onto the nuclei structure
        for (int jnuc=0; jnuc<=recvnonucs-1; jnuc++){
            nucleus* temp = new nucleus;
            temp->index=r_index[jnuc];
            dealii::Point<problemDIM> r_center;
            r_center[0]=r_center_x[jnuc];
            r_center[1]=r_center_y[jnuc];
            temp->center=r_center;
            temp->radius=r_radius[jnuc];
            temp->seededTime=r_seededTime[jnuc];
            temp->seedingTime = r_seedingTime[jnuc];
            temp->seedingTimestep = r_seedingTimestep[jnuc];
            newnuclei.push_back(*temp);
        }

    }
}

void broadcastUpdate (std::vector<nucleus> &newnuclei, int broadcastProc, int thisProc)
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
        //4 - radius
        std::vector<double> r_radius(initial_vec_size,0.0);
        //5 - seededTime
        std::vector<double> r_seededTime(initial_vec_size,0.0);
        //6 - seedingTime
        std::vector<double> r_seedingTime(initial_vec_size,0.0);
        //7 - seedingTimestep
        std::vector<unsigned int> r_seedingTimestep(initial_vec_size,0);
        
        if (thisProc == broadcastProc){
        	for (std::vector<nucleus>::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){
        		r_index.push_back(thisNuclei->index);
        		dealii::Point<problemDIM> s_center=thisNuclei->center;
        		r_center_x.push_back(s_center[0]);
        		r_center_y.push_back(s_center[1]);
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
            temp->center=r_center;
            temp->radius=r_radius[jnuc];
            temp->seededTime=r_seededTime[jnuc];
            temp->seedingTime = r_seedingTime[jnuc];
            temp->seedingTimestep = r_seedingTimestep[jnuc];
            //nuclei.push_back(*temp); // I want to save this for later -Steve
            newnuclei.push_back(*temp);
        }
        
    }
}

void resolveNucleationConflicts (std::vector<nucleus> &newnuclei){

	std::vector<nucleus> newnuclei_cleaned;
	unsigned int old_num_nuclei = nuclei.size();

	for (unsigned int nuc_index=0; nuc_index<newnuclei.size(); nuc_index++){
		bool isClose=false;

		for (unsigned int prev_nuc_index=0; prev_nuc_index<nuc_index; prev_nuc_index++){

			// We may want to break this section into a separate function to allow different choices for when
			// nucleation should be prevented
			if (newnuclei[nuc_index].center.distance(newnuclei[prev_nuc_index].center) < minDistBetwenNuclei){
				isClose = true;
				std::cout << "Conflict between nuclei! Distance is: " << newnuclei[nuc_index].center.distance(newnuclei[prev_nuc_index].center) << std::endl;
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
// Global nucleation procedure
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::modifySolutionFields()
{
	//current time step
	unsigned int inc=this->currentIncrement;

    if ( (inc <= timeIncrements) && (inc % skipNucleationSteps == 0) ){
    	//current time
		double t=this->currentTime;

		//get the list of node points in the domain
		std::map<dealii::types::global_dof_index, dealii::Point<dim> > support_points;
		dealii::DoFTools::map_dofs_to_support_points (dealii::MappingQ1<dim>(), *this->dofHandlersSet[0], support_points);
		//fields
		vectorType* n=this->solutionSet[this->getFieldIndex("n")];
		vectorType* c=this->solutionSet[this->getFieldIndex("c")];


        //vector of all the NEW nuclei seeded in this time step
        std::vector<nucleus> newnuclei;
        
        // Attempt to nucleate (each processor independently)
        nucAttempt(newnuclei, support_points, c, n, t, inc);

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

        // Add the new nuclei to the list of nuclei
        nuclei.insert(nuclei.end(),newnuclei.begin(),newnuclei.end());

    }
}
