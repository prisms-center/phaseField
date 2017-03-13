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
	double J=k1*exp(-k2/(std::max(cValue-calmin,1.0e-6)));
	double retProb=1.0-exp(-J*timeStep*((double)skipNucleationSteps)*dV);
    return retProb;
}

// =================================================================================
// Get list of prospective new nuclei for the local processor
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::getLocalNucleiList(std::vector<nucleus<dim>> &newnuclei) const
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
                bool insafetyzone = true;
                for (unsigned int j=0; j < dim; j++){
                    std::string BC_type = this->BC_list[1].var_BC_type[2*j];
                    bool periodic_j = (BC_type=="PERIODIC");
                    bool insafetyzone_j = (periodic_j || ((q_point_list[q_point][j] > borderreg) && (q_point_list[q_point][j] < this->userInputs.domain_size[j]-borderreg)));
                	insafetyzone = insafetyzone && insafetyzone_j;
                }
				if (insafetyzone){
					if (var_value2[q_point] < maxOrderParameterNucleation){

						//Compute random no. between 0 and 1 (old method)
						rand_val=distr(gen);
						double Prob=nucProb(var_value[q_point],fe_values.JxW(q_point));
						if (rand_val <= Prob){
							//loop over all existing nuclei to check if they are in the vicinity
							bool isClose=false;

							// ---------------------------------------------------------------------------------
							// Check to see if the prospective nucleus would overlap with any other particles
				    			typename DoFHandler<dim>::active_cell_iterator   di_overlap = this->dofHandlersSet_nonconst[0]->begin_active();
				    			while (di_overlap != this->dofHandlersSet_nonconst[0]->end())
				    			{
				    				if (di_overlap->is_locally_owned()){
				    					fe_values.reinit (di_overlap);
				    					fe_values.get_function_values(*(this->solutionSet[1]), var_value2_overlap);
				    					q_point_list_overlap = fe_values.get_quadrature_points();
				    					for (unsigned int q_point_overlap=0; q_point_overlap<num_quad_points; ++q_point_overlap){
				    						double weighted_dist = 0.0;
				    						for (unsigned int i=0; i<dim; i++){
				    							double temp = (q_point_list[q_point](i)-q_point_list_overlap[q_point_overlap](i))/opfreeze_semiaxes[i];
				    							weighted_dist += temp*temp;
				    						}
				    						if (weighted_dist < 1.0){
				    							if (var_value2_overlap[q_point_overlap] > 0.1){
				    								isClose=true;
				    								std::cout << "Attempted nucleation failed due to overlap w/ existing particle!!!!!!" << q_point_list_overlap[q_point_overlap](0) << " " << q_point_list_overlap[q_point_overlap](1) << " " << q_point_list[q_point](0) << " " << q_point_list[q_point](1) << std::endl;
				    								break;
				    							}
				    						}
				    					}
				    					if (isClose) break;
				    				}
				    				++di_overlap;
				    			}
				    			fe_values.reinit (di);
							// ---------------------------------------------------------------------------------

							if (!isClose){
								std::cout << "Nucleation event. Nucleus no. " << nuclei.size()+1 << std::endl;
								std::cout << "nucleus center " << q_point_list[q_point] << std::endl;
								nucleus<dim>* temp = new nucleus<dim>;
								temp->index=nuclei.size();
								temp->center=q_point_list[q_point];
								temp->semiaxes.push_back(semiaxis_a);
								temp->semiaxes.push_back(semiaxis_b);
								if (dim == 3){
									temp->semiaxes.push_back(semiaxis_c);
								}
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
// Refine mesh near the new nuclei
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::refineMeshNearNuclei(std::vector<nucleus<dim>> newnuclei)
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
					for (typename std::vector<nucleus<dim>>::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){
						double weighted_dist = 0.0;
						for (unsigned int i=0; i<dim; i++){
							double temp = (thisNuclei->center(i) - q_point_list[q_point](i))/(opfreeze_semiaxes[i]);
							weighted_dist += temp*temp;
						}
						if (weighted_dist < 1.0){
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
		std::vector<nucleus<dim>> newnuclei;

		// Get list of prospective new nuclei for the local processor
		this->pcout << "Nucleation attempt for increment " << this->currentIncrement << std::endl;
		getLocalNucleiList(newnuclei);

		// Generate global list of new nuclei and resolve conflicts between new nuclei
		parallelNucleationList<dim> new_nuclei_parallel(newnuclei);
		newnuclei = new_nuclei_parallel.buildGlobalNucleiList(minDistBetweenNuclei, nuclei.size());

        // Add the new nuclei to the list of nuclei
        nuclei.insert(nuclei.end(),newnuclei.begin(),newnuclei.end());

        // Refine mesh near the new nuclei
        if (newnuclei.size() > 0 && this->userInputs.h_adaptivity == true){
        	refineMeshNearNuclei(newnuclei);
        }
    }
}
