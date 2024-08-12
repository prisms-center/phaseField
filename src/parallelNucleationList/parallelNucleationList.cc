#include "../../include/parallelNucleationList.h"
#include "../../include/nucleus.h"
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>
#include <iostream>

// =================================================================================
// Constructor
// =================================================================================
template <int dim>
parallelNucleationList<dim>::parallelNucleationList (std::vector<nucleus<dim> > _newnuclei): newnuclei(_newnuclei) {}

// =================================================================================
// Generate global list of new nuclei and resolve conflicts between new nuclei
// =================================================================================
template <int dim>
std::vector<nucleus<dim> > parallelNucleationList<dim>::buildGlobalNucleiList(double min_dist_between_nuclei, double min_dist_between_OP, unsigned int old_num_nuclei)
{
	//MPI INITIALIZATON
	int numProcs=dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
	int thisProc=dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
	if (numProcs > 1) {
		// Cycle through each processor, sending and receiving, to append the list of new nuclei
		for (int proc_index=0; proc_index < numProcs-1; proc_index++){
			if (thisProc == proc_index){
				sendUpdate(thisProc+1);
			}
			else if (thisProc == proc_index+1){
				receiveUpdate(thisProc-1);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		// The final processor now has all of the new nucleation attempts
		// Check for conflicts on the final processor before broadcasting the list
		if (thisProc == numProcs-1){
			resolveNucleationConflicts(min_dist_between_nuclei, min_dist_between_OP, old_num_nuclei);
		}

		// The final processor now has the final list of the new nuclei, broadcast it to all the other processors
		broadcastUpdate(numProcs-1, thisProc);
	}
	else {
		// Check for conflicts between nucleation attempts this time step
		resolveNucleationConflicts(min_dist_between_nuclei, min_dist_between_OP, old_num_nuclei);
	}

	return newnuclei;
}

// =================================================================================
// Sends the list of new nuclei to the next processor
// =================================================================================
template <int dim>
void parallelNucleationList<dim>::sendUpdate (int procno) const
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
        std::vector<double> s_semiaxis_a;
        std::vector<double> s_semiaxis_b;
        std::vector<double> s_semiaxis_c;
        //6 - seededTime
        std::vector<double> s_seededTime;
        //7 - seedingTime
        std::vector<double> s_seedingTime;
        //8 - seedingTimestep
        std::vector<unsigned int> s_seedingTimestep;
        //9 - orderParameterIndex
        std::vector<unsigned int> s_orderParameterIndex;

        //Loop to store info of all nuclei into vectors
        for (typename std::vector<nucleus<dim> >::const_iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){
            s_index.push_back(thisNuclei->index);
            dealii::Point<dim> s_center=thisNuclei->center;
            s_center_x.push_back(s_center[0]);
            s_center_y.push_back(s_center[1]);
            if (dim ==3)
            	s_center_z.push_back(s_center[2]);

            s_semiaxis_a.push_back(thisNuclei->semiaxes[0]);
            s_semiaxis_b.push_back(thisNuclei->semiaxes[1]);
            if (dim == 3)
            	s_semiaxis_c.push_back(thisNuclei->semiaxes[2]);

            s_seededTime.push_back(thisNuclei->seededTime);
            s_seedingTime.push_back(thisNuclei->seedingTime);
            s_seedingTimestep.push_back(thisNuclei->seedingTimestep);
            s_orderParameterIndex.push_back(thisNuclei->orderParameterIndex);
        }
        //Send vectors to next processor
        MPI_Send(&s_index[0], currnonucs, MPI_UNSIGNED, procno, 1, MPI_COMM_WORLD);
        MPI_Send(&s_center_x[0], currnonucs, MPI_DOUBLE, procno, 2, MPI_COMM_WORLD);
        MPI_Send(&s_center_y[0], currnonucs, MPI_DOUBLE, procno, 3, MPI_COMM_WORLD);
        if (dim ==3){
        	MPI_Send(&s_center_z[0], currnonucs, MPI_DOUBLE, procno, 4, MPI_COMM_WORLD);
        	MPI_Send(&s_semiaxis_c[0], currnonucs, MPI_DOUBLE, procno, 7, MPI_COMM_WORLD);
        }

        MPI_Send(&s_semiaxis_a[0], currnonucs, MPI_DOUBLE, procno, 5, MPI_COMM_WORLD);
        MPI_Send(&s_semiaxis_b[0], currnonucs, MPI_DOUBLE, procno, 6, MPI_COMM_WORLD);

        MPI_Send(&s_seededTime[0], currnonucs, MPI_DOUBLE, procno, 8, MPI_COMM_WORLD);
        MPI_Send(&s_seedingTime[0], currnonucs, MPI_DOUBLE, procno, 9, MPI_COMM_WORLD);
        MPI_Send(&s_seedingTimestep[0], currnonucs, MPI_UNSIGNED, procno, 10, MPI_COMM_WORLD);
        MPI_Send(&s_orderParameterIndex[0], currnonucs, MPI_UNSIGNED, procno, 11, MPI_COMM_WORLD);
    }
    //END OF MPI SECTION
}

// =================================================================================
// Recieves the list of new nuclei to the next processor
// =================================================================================
template <int dim>
void parallelNucleationList<dim>::receiveUpdate (int procno)
{
    //MPI PROCEDURE TO RECIEVE INFORMATION FROM ANOTHER PROCESSOR AND UPDATE LOCAL NUCLEI INFORMATION
    int recvnonucs = 0;
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
        //5 - semiaxes
        std::vector<double> r_semiaxis_a(recvnonucs,0.0);
        std::vector<double> r_semiaxis_b(recvnonucs,0.0);
        std::vector<double> r_semiaxis_c(recvnonucs,0.0);
        //6 - seededTime
        std::vector<double> r_seededTime(recvnonucs,0.0);
        //7 - seedingTime
        std::vector<double> r_seedingTime(recvnonucs,0.0);
        //8 - seedingTimestep
        std::vector<unsigned int> r_seedingTimestep(recvnonucs,0);
        //9 - orderParameterIndex
        std::vector<unsigned int> r_orderParameterIndex(recvnonucs,0);

        //Recieve vectors from processor procno
        MPI_Recv(&r_index[0], recvnonucs, MPI_UNSIGNED, procno, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_center_x[0], recvnonucs, MPI_DOUBLE, procno, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_center_y[0], recvnonucs, MPI_DOUBLE, procno, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (dim ==3){
        	MPI_Recv(&r_center_z[0], recvnonucs, MPI_DOUBLE, procno, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	MPI_Recv(&r_semiaxis_c[0], recvnonucs, MPI_DOUBLE, procno, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Recv(&r_semiaxis_a[0], recvnonucs, MPI_DOUBLE, procno, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_semiaxis_b[0], recvnonucs, MPI_DOUBLE, procno, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Recv(&r_seededTime[0], recvnonucs, MPI_DOUBLE, procno, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_seedingTime[0], recvnonucs, MPI_DOUBLE, procno, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_seedingTimestep[0], recvnonucs, MPI_UNSIGNED, procno, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&r_orderParameterIndex[0], recvnonucs, MPI_UNSIGNED, procno, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //Loop to store info in vectors onto the nuclei structure
        for (int jnuc=0; jnuc<=recvnonucs-1; jnuc++){
            nucleus<dim>* temp = new nucleus<dim>;
            temp->index=r_index[jnuc];
            dealii::Point<dim> r_center;
            r_center[0]=r_center_x[jnuc];
            r_center[1]=r_center_y[jnuc];
            if (dim ==3)
            	r_center[2]=r_center_z[jnuc];
            temp->center=r_center;
            temp->semiaxes.push_back(r_semiaxis_a[jnuc]);
            temp->semiaxes.push_back(r_semiaxis_b[jnuc]);
            if (dim ==3)
            	temp->semiaxes.push_back(r_semiaxis_c[jnuc]);
            temp->seededTime=r_seededTime[jnuc];
            temp->seedingTime = r_seedingTime[jnuc];
            temp->seedingTimestep = r_seedingTimestep[jnuc];
            temp->orderParameterIndex = r_orderParameterIndex[jnuc];
            newnuclei.push_back(*temp);
        }

    }
}

// =================================================================================
// Broadcast the final list of new nuclei from the last processor to the rest
// =================================================================================
template <int dim>
void parallelNucleationList<dim>::broadcastUpdate (int broadcastProc, int thisProc)
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
        std::vector<double> r_semiaxis_a(initial_vec_size,0.0);
        std::vector<double> r_semiaxis_b(initial_vec_size,0.0);
        std::vector<double> r_semiaxis_c(initial_vec_size,0.0);
        //6 - seededTime
        std::vector<double> r_seededTime(initial_vec_size,0.0);
        //7 - seedingTime
        std::vector<double> r_seedingTime(initial_vec_size,0.0);
        //8 - seedingTimestep
        std::vector<unsigned int> r_seedingTimestep(initial_vec_size,0);
        //9 - orderParameterIndex
        std::vector<unsigned int> r_orderParameterIndex(initial_vec_size,0);

        if (thisProc == broadcastProc){
        	for (typename std::vector<nucleus<dim> >::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){
        		r_index.push_back(thisNuclei->index);
        		dealii::Point<dim> s_center=thisNuclei->center;
        		r_center_x.push_back(s_center[0]);
        		r_center_y.push_back(s_center[1]);
                if (dim ==3){
                	r_center_z.push_back(s_center[2]);
				}

                r_semiaxis_a.push_back(thisNuclei->semiaxes[0]);
                r_semiaxis_b.push_back(thisNuclei->semiaxes[1]);
                if (dim ==3){
                	r_semiaxis_c.push_back(thisNuclei->semiaxes[2]);
				}

        		r_seededTime.push_back(thisNuclei->seededTime);
        		r_seedingTime.push_back(thisNuclei->seedingTime);
        		r_seedingTimestep.push_back(thisNuclei->seedingTimestep);
        		r_orderParameterIndex.push_back(thisNuclei->orderParameterIndex);
        	}
        }

        //Recieve vectors from processor procno
        MPI_Bcast(&r_index[0], currnonucs, MPI_UNSIGNED, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_center_x[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_center_y[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        if (dim ==3)
        	MPI_Bcast(&r_center_z[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);

        MPI_Bcast(&r_semiaxis_a[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_semiaxis_b[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        if (dim == 3)
        	MPI_Bcast(&r_semiaxis_c[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);

        MPI_Bcast(&r_seededTime[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_seedingTime[0], currnonucs, MPI_DOUBLE, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_seedingTimestep[0], currnonucs, MPI_UNSIGNED, broadcastProc, MPI_COMM_WORLD);
        MPI_Bcast(&r_orderParameterIndex[0], currnonucs, MPI_UNSIGNED, broadcastProc, MPI_COMM_WORLD);

        newnuclei.clear();

        //Loop to store info in vectors onto the nuclei structure
        for (int jnuc=0; jnuc<=currnonucs-1; jnuc++){
            nucleus<dim>* temp = new nucleus<dim>;
            temp->index=r_index[jnuc];
            dealii::Point<dim> r_center;
            r_center[0]=r_center_x[jnuc];
            r_center[1]=r_center_y[jnuc];
            if (dim ==3)
            	r_center[2]=r_center_z[jnuc];
            temp->center=r_center;
            temp->semiaxes.push_back(r_semiaxis_a[jnuc]);
            temp->semiaxes.push_back(r_semiaxis_b[jnuc]);
            if (dim == 3){
            	temp->semiaxes.push_back(r_semiaxis_c[jnuc]);
            }
            temp->seededTime=r_seededTime[jnuc];
            temp->seedingTime = r_seedingTime[jnuc];
            temp->seedingTimestep = r_seedingTimestep[jnuc];
            temp->orderParameterIndex = r_orderParameterIndex[jnuc];
            newnuclei.push_back(*temp);
        }

    }
}

// =================================================================================
// Determine if any new nuclei are in conflict and resolve those conflicts
// =================================================================================
template <int dim>
void parallelNucleationList<dim>::resolveNucleationConflicts (double min_dist_between_nuclei, double min_dist_between_OP, unsigned int old_num_nuclei)
{
	std::vector<nucleus<dim> > newnuclei_cleaned;

	for (unsigned int nuc_index=0; nuc_index<newnuclei.size(); nuc_index++){
		bool isClose=false;

		for (unsigned int prev_nuc_index=0; prev_nuc_index<nuc_index; prev_nuc_index++){

			// We may want to break this section into a separate function to allow different choices for when
			// nucleation should be prevented
		
		// This section checks whether there is any overlap between the two nuclei , irrespective of there order parameter
	    if (newnuclei[nuc_index].center.distance(newnuclei[prev_nuc_index].center) < min_dist_between_nuclei){
				isClose = true;
				std::cout << "Conflict between nuclei! Distance is: " << newnuclei[nuc_index].center.distance(newnuclei[prev_nuc_index].center)
						<< " Conflict removed."<< std::endl;
				break;
			} else if(newnuclei[nuc_index].orderParameterIndex == newnuclei[prev_nuc_index].orderParameterIndex && newnuclei[nuc_index].center.distance(newnuclei[prev_nuc_index].center) < min_dist_between_OP){ // This section makes sure two nuclei with same order parameter are located further apart to avoid any coalesence during growth

			    	isClose = true;
				std::cout << "Conflict between nuclei! Distance is: " << newnuclei[nuc_index].center.distance(newnuclei[prev_nuc_index].center)
						<< " Conflict removed."<< std::endl;
				break;
			}
		}

		if (!isClose){
			newnuclei[nuc_index].index = old_num_nuclei + newnuclei_cleaned.size();
			std::cout << "Nuclei number : " << (old_num_nuclei + newnuclei_cleaned.size()) << std::endl;
			newnuclei_cleaned.push_back(newnuclei[nuc_index]);
		}
	}

	newnuclei = newnuclei_cleaned;
}

// =================================================================================
// Remove nuclei from the list of nuclei given a local list of nucleus indices
// =================================================================================
template <int dim>
std::vector<nucleus<dim> > parallelNucleationList<dim>::removeSubsetOfNuclei(std::vector<unsigned int> nuclei_to_remove, unsigned int nuclei_size){
	// Note: This method is very similar to buildGlobalNucleiList in structure, and uses simplified versions of what is done
	// in sendUpdate, receiveUpdate, and broadcastUpdate. There is likely a cleaner way to reorganize the methods to reduce
	// duplication.

	//MPI INITIALIZATON
	int numProcs=dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
	int thisProc=dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

	// Build a global list of nuclei to delete, first sending the length of the vector of indices, then the vector itself
	if (numProcs > 1) {
		// Cycle through each processor, sending and receiving, to append the list of new nuclei
		for (int proc_index=0; proc_index < numProcs-1; proc_index++){
			if (thisProc == proc_index){
				int currnonucs= nuclei_to_remove.size();
				MPI_Send(&currnonucs, 1, MPI_INT, thisProc+1, 0, MPI_COMM_WORLD);
				MPI_Send(&nuclei_to_remove[0], currnonucs, MPI_UNSIGNED, thisProc+1, 1, MPI_COMM_WORLD);
			}
			else if (thisProc == proc_index+1){
				int recvnonucs = 0;
				MPI_Recv(&recvnonucs, 1, MPI_INT, thisProc-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				std::vector<unsigned int> recieved_nuclei_to_remove(recvnonucs);
				MPI_Recv(&recieved_nuclei_to_remove[0], recvnonucs, MPI_UNSIGNED, thisProc-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				nuclei_to_remove.insert(nuclei_to_remove.end(), recieved_nuclei_to_remove.begin(), recieved_nuclei_to_remove.end());

			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		// The final processor now has the final list of the new nuclei, broadcast it to all the other processors
		int currnonucs= nuclei_to_remove.size();
		MPI_Bcast(&currnonucs, 1, MPI_INT, numProcs-1, MPI_COMM_WORLD);
		std::vector<unsigned int> recieved_nuclei_to_remove(currnonucs);
		if (thisProc == numProcs-1){
			recieved_nuclei_to_remove = nuclei_to_remove;
		}
		MPI_Bcast(&recieved_nuclei_to_remove[0], currnonucs, MPI_UNSIGNED, numProcs-1, MPI_COMM_WORLD);
		nuclei_to_remove = recieved_nuclei_to_remove;
	}

	for (unsigned int i=0; i<nuclei_to_remove.size(); i++){
		std::cout << thisProc << ": " << nuclei_to_remove[i] << std::endl;
	}

	// Remove the nuclei from the list
	std::vector<nucleus<dim> > pruned_list;
	for (unsigned int nuc = 0; nuc <  newnuclei.size(); nuc++){
		bool pruneNucleus = false;
		for (unsigned int i = 0; i < nuclei_to_remove.size(); i++){
			if (nuclei_to_remove[i] == nuclei_size+nuc){
				pruneNucleus = true;
				break;
			}
		}
		if (!pruneNucleus){
			pruned_list.push_back(newnuclei[nuc]);
		}
	}
	return pruned_list;
}

// =================================================================================
// Template instantiations
// =================================================================================
template class parallelNucleationList<2>;
template class parallelNucleationList<3>;
