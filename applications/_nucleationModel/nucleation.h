// =================================================================================
//Function prototypes
// =================================================================================
double nucProb(double cValue);
void nucAttempt(std::vector<nucleus> &newnuclei, std::map<dealii::types::global_dof_index, dealii::Point<problemDIM> > support_points, vectorType* c, vectorType* n, double t, unsigned int inc);
void receiveUpdate (std::vector<nucleus> &newnuclei, int procno);
void sendUpdate (int procno);

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
	double J=k1*exp(-k2/(cValue-calmin));
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
                        nuclei.push_back(*temp);
                        newnuclei.push_back(*temp);
                    }
                }
            }
        }
    }
}

void sendUpdate (int procno)
{
    int currnonucs=nuclei.size();
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
        for (std::vector<nucleus>::iterator thisNuclei=nuclei.begin(); thisNuclei!=nuclei.end(); ++thisNuclei){
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
        MPI_Send(&s_index[0], currnonucs, MPI_UNSIGNED, procno+1, 1, MPI_COMM_WORLD);
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
    int currnonucs = nuclei.size();
    MPI_Recv(&recvnonucs, 1, MPI_INT, procno, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (recvnonucs > currnonucs){
        
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
        for (int jnuc=currnonucs; jnuc<=recvnonucs-1; jnuc++){
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
            nuclei.push_back(*temp);
            newnuclei.push_back(*temp);
        }
        
    }
}

// =================================================================================
// Global nucleation procedure
// =================================================================================
template <int dim>
void generalizedProblem<dim>::modifySolutionFields()
{
    //current time
    double t=this->currentTime;
    //current time step
    unsigned int inc=this->currentIncrement;
    
    //get the list of node points in the domain
    std::map<dealii::types::global_dof_index, dealii::Point<dim> > support_points;
    dealii::DoFTools::map_dofs_to_support_points (dealii::MappingQ1<dim>(), *this->dofHandlersSet[0], support_points);
    //fields
    vectorType* n=this->solutionSet[this->getFieldIndex("n")];
    vectorType* c=this->solutionSet[this->getFieldIndex("c")];
    
    if ( (inc <= timeIncrements) && (inc % skipNucleationSteps == 0) ){
        
        //vector of all the NEW nuclei seeded in this time step
        std::vector<nucleus> newnuclei;
        
        //MPI INITIALIZATON
        int numProcs=Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
        int thisProc=Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
        
        if (numProcs == 1) {
        	//Serial option
            nucAttempt(newnuclei, support_points, c, n, t, inc);
        } else {
        	//Parallel option
            //All processors except 0 receive updated info from previous processor
            if (thisProc > 0){
                //MPI SECTION TO RECIEVE INFORMATION FROM THE PREVIOUS PROCESSOR
                receiveUpdate(newnuclei, thisProc-1);
                std::cout << "updated info recieved by next proc." << std::endl;
            }
            //Each processor attempts nucleation
            if (thisProc == 0){
            	//PLEASE UNCOMMENT THIS LINE WHEN CODE IS FIXED (DM)
                //std::cout << "nucleation attempt" << std::endl;
            }
            nucAttempt(newnuclei, support_points, c, n, t, inc);
            //All processors (except last) send info to next processor
            if (thisProc < numProcs-1){
                //MPI SECTION TO SEND INFORMATION TO THE NEXT PROCESSOR
                sendUpdate(thisProc+1);
                std::cout << "updated info sent to next proc." << std::endl;
            }
            //Last processor (N-1) sends info to processor 0
            if (thisProc == numProcs-1){
                //MPI SECTION TO SEND INFORMATION TO PROCESSOR 0
                sendUpdate(0);
                std::cout << "updated info sent to proc. 0" << std::endl;
            }
            if (thisProc == 0){
                //MPI SECTION TO RECIEVE INFORMATION FROM PROCESSOR N-1
                receiveUpdate(newnuclei, numProcs-1);
                std::cout << "updated info recieved by to proc. 0" << std::endl;
            }
            //Barrier for all processors
            MPI_Barrier(MPI_COMM_WORLD);
            std::cout << "Passed 1st MPI barrier" << std::endl;
            
            if (thisProc == 0){
                //MPI SECTION TO BROADCAST INFORMATION FROM PROCESSOR 0 TO ALL OTHER PROCESSORS
                //Sending local no. of nuclei
                for (int jproc=1; jproc<=numProcs-1; jproc++){
                    sendUpdate(jproc);
                }
                std::cout << "Broadcast from 0 to all procs" << std::endl;
            }
            if (thisProc > 0){
                //MPI SECTION TO RECIEVE INFORMATION FROM THE PROCESSOR 0
                receiveUpdate(newnuclei, 0);
                std::cout << "Local proc. recieves broadcast from 0" << std::endl;
            }
            //Barrier for all processor
            MPI_Barrier(MPI_COMM_WORLD);
            std::cout << "Passed 2nd MPI barrier" << std::endl;
        }
        //Seeding nucleus section
        //Looping over all nodes
        for (typename std::map<dealii::types::global_dof_index, dealii::Point<dim> >::iterator it=support_points.begin(); it!=support_points.end(); ++it){
            unsigned int dof=it->first;
            //set only local owned values of the parallel vector
            if (n->locally_owned_elements().is_element(dof)){
                dealii::Point<dim> nodePoint=it->second;
                //Looping over all nuclei seeded in this iteration
                for (std::vector<nucleus>::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){
                    dealii::Point<dim> center=thisNuclei->center;
                    double r=nodePoint.distance(center);
                    if (r<=opfreeze_radius){
                        (*n)(dof)=0.5*(1.0-std::tanh((r-n_radius)/interface_coeff));
                    }
                }
            }
        }
    }
}
