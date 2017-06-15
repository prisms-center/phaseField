// =================================================================================
// NUCLEATION FUNCTIONS
// =================================================================================

// =================================================================================
// Nucleation probability
// =================================================================================
template <int dim, int degree>
double customPDE<dim,degree>::nucProb(double cValue, double dV, double ct) const
{
	//Supersaturation factor
    double ssf;
    if (dim ==2) ssf=cValue-calmin;
    if (dim ==3) ssf=(cValue-calmin)*(cValue-calmin);
	// Calculate the nucleation rate
	double J=k1*exp(-k2/(std::max(ssf,1.0e-6)))*exp(-tau/ct);
	double retProb=1.0-exp(-J*userInputs.dtValue*((double)skipNucleationSteps)*dV);
    return retProb;
}

// =================================================================================
// Get list of prospective new nuclei for the local processor
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::getLocalNucleiList(std::vector<nucleus<dim> > &newnuclei) const
{
	// Nickname for current time and time step
	double t=this->currentTime;
	unsigned int inc=this->currentIncrement;

    //QGauss<dim>  quadrature(degree+1);
    QGaussLobatto<dim>  quadrature(degree+1);
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

	//Element cycle
	while (di != this->dofHandlersSet_nonconst[0]->end())
	{
		if (di->is_locally_owned()){
            //Obtaining average element concentration by averaging over element's quadrature points
            fe_values.reinit(di);
            fe_values.get_function_values(*(this->solutionSet[0]), var_value);
            fe_values.get_function_values(*(this->solutionSet[1]), var_value2);
            q_point_list = fe_values.get_quadrature_points();
            double ele_vol = 0.0;
            double ele_av_conc = 0.0;
        	// Loop over the quadrature points
            for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                ele_av_conc = ele_av_conc + var_value[q_point]*fe_values.JxW(q_point);
                ele_vol = ele_vol + fe_values.JxW(q_point);
            }
            ele_av_conc = ele_av_conc/ele_vol;

            //Compute random no. between 0 and 1 (new method)
            rand_val=distr(gen);
            //Nucleation probability
            double Prob=nucProb(ele_av_conc,ele_vol,t);

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
                    bool insafetyzone_j = (periodic_j || ((nuc_ele_pos[j] > borderreg) && (nuc_ele_pos[j] < userInputs.domain_size[j]-borderreg)));
                    insafetyzone = insafetyzone && insafetyzone_j;
                }

                if (insafetyzone){

                	// Check to see if the order parameter anywhere within the element is above the threshold
                    bool anyqp_OK = false;
                    for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                        if (var_value2[q_point] < maxOrderParameterNucleation){
                            anyqp_OK =true;
                        }
                    }

                    if (anyqp_OK){
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
void customPDE<dim,degree>::safetyCheckNewNuclei(std::vector<nucleus<dim> > newnuclei, std::vector<unsigned int> &conflict_inds)
{
    //QGauss<dim>  quadrature(degree+1);
    QGaussLobatto<dim>  quadrature(degree+1);
    FEValues<dim> fe_values (*(this->FESet[0]), quadrature, update_values|update_quadrature_points|update_JxW_values);
    const unsigned int   num_quad_points = quadrature.size();
    std::vector<double> var_value2(num_quad_points);
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
                fe_values.get_function_values(*(this->solutionSet[1]), var_value2);
                q_point_list = fe_values.get_quadrature_points();

                //Quadrature points cycle
                for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                    // Calculate the ellipsoidal distance to the center of the nucleus
                    double weighted_dist = 0.0;
                    for (unsigned int i=0; i<dim; i++){
                        double shortest_edist = thisNuclei->center(i) - q_point_list[q_point](i);
                        bool periodic_i = (userInputs.BC_list[1].var_BC_type[2*i]==PERIODIC);
                        if (periodic_i){
                            double domsize =userInputs.domain_size[i];
                            shortest_edist = shortest_edist-round(shortest_edist/domsize)*domsize;
                        }
                        double temp = shortest_edist/(opfreeze_semiaxes[i]);
                        weighted_dist += temp*temp;
                    }
                    if (weighted_dist < 1.0){
                        if (var_value2[q_point] > 0.1){
                            isClose=true;
                            std::cout << "Attempted nucleation failed due to overlap w/ existing particle!!!!!!"  << std::endl;
                            std::cout << "Nucleus index " << thisNuclei->index << std::endl;
                            std::cout << "Nucleus center " << thisNuclei->center << std::endl;
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

    for (unsigned int remesh_index=0; remesh_index < (userInputs.max_refinement_level-userInputs.min_refinement_level); remesh_index++){
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
                    diag_dist += (userInputs.domain_size[i]*userInputs.domain_size[i])/(userInputs.subdivisions[i]*userInputs.subdivisions[i]);
                }
                diag_dist = sqrt(diag_dist);
                diag_dist /= 2.0*pow(2.0,ti->level());

                for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
                    for (typename std::vector<nucleus<dim> >::iterator thisNuclei=newnuclei.begin(); thisNuclei!=newnuclei.end(); ++thisNuclei){

                        // Calculate the ellipsoidal distance to the center of the nucleus
                        double weighted_dist = 0.0;
                        for (unsigned int i=0; i<dim; i++){
                            double shortest_edist = thisNuclei->center(i) - q_point_list[q_point](i);
                            bool periodic_i = (userInputs.BC_list[1].var_BC_type[2*i]==PERIODIC);
                            if (periodic_i){
                                double domsize =userInputs.domain_size[i];
                                shortest_edist = shortest_edist-round(shortest_edist/domsize)*domsize;
                            }
                            double temp = shortest_edist/(opfreeze_semiaxes[i]);
                            weighted_dist += temp*temp;
                        }

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
        std::vector<nucleus<dim> > newnuclei;

        // Get list of prospective new nuclei for the local processor
        this->pcout << "Nucleation attempt for increment " << this->currentIncrement << std::endl;
        getLocalNucleiList(newnuclei);

        // Generate global list of new nuclei and resolve conflicts between new nuclei
        parallelNucleationList<dim> new_nuclei_parallel(newnuclei);
        newnuclei = new_nuclei_parallel.buildGlobalNucleiList(minDistBetweenNuclei, nuclei.size());

        // Final check to resolve overlap conflicts with existing precipitates
        std::vector<unsigned int> conflict_inds;
        safetyCheckNewNuclei(newnuclei, conflict_inds);

        newnuclei = new_nuclei_parallel.removeSubsetOfNuclei(conflict_inds, nuclei.size());

        // Add the new nuclei to the list of nuclei
        nuclei.insert(nuclei.end(),newnuclei.begin(),newnuclei.end());

        // Refine mesh near the new nuclei
        if (newnuclei.size() > 0 && userInputs.h_adaptivity == true){
            refineMeshNearNuclei(newnuclei);
        }
        //Print total no. of nuclei after nucleation attempt
        this->pcout << "Print total no. of nuclei after nucleation attempt" << std::endl;
        this->pcout << "Increment: " << this->currentIncrement << " No. nuclei: " <<  nuclei.size() << std::endl;
    }
}
