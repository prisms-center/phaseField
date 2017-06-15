// Methods in MatrixFreePDE to update the list of nuclei

template <int dim, int degree>
void MatrixFreePDE<dim,degree>::updateNucleiList(){

    if (userInputs.nucleation_occurs){
        if ( this->currentIncrement % skipNucleationSteps == 0 ){
            computing_timer.enter_section("matrixFreePDE: nucleation");

            if (this->currentIncrement == 1){
                while (newnuclei.size() == 0){
                    this->currentTime+=this->userInputs.dtValue*(double)skipNucleationSteps;
                    this->currentIncrement+=skipNucleationSteps;

                    while (this->outputTimeStepList.size() > 0 && this->outputTimeStepList[0] < this->currentIncrement){
                        this->outputTimeStepList.erase(this->outputTimeStepList.begin());
                    }

                    std::vector<nucleus> new_nuclei = getNewNuclei();
                }
            }
            else {
                std::vector<nucleus> new_nuclei = getNewNuclei();
            }

            nuclei.insert(nuclei.end(),new_nuclei.begin(),new_nuclei.end());

            if (new_nuclei.size() > 0 && userInputs.h_adaptivity == true){
                refineMeshNearNuclei(newnuclei);
            }

            computing_timer.exit_section("matrixFreePDE: nucleation");
        }
    }

}

template <int dim, int degree>
std::vector<nucleus> MatrixFreePDE<dim,degree>::getNewNuclei(){

    // Declare a vector of all the NEW nuclei seeded in this time step
    std::vector<nucleus<dim> > newnuclei;

    // Get list of prospective new nuclei for the local processor
    this->pcout << "Nucleation attempt for increment " << this->currentIncrement << std::endl;

    // Move something like these to userInputParameters
    // std::vector<unsigned int> order_parameter_list;
    // std::vector<unsigned int> other_var_list;
    // other_var_list.push_back(0);
    // order_parameter_list.push_back(1);
    // order_parameter_list.push_back(2);
    // order_parameter_list.push_back(3);

    getLocalNucleiList(newnuclei,order_parameter_list,other_var_list);
    this->pcout << "nucleation attempt! " << this->currentTime << " " << this->currentIncrement << std::endl;

    // Generate global list of new nuclei and resolve conflicts between new nuclei
    parallelNucleationList<dim> new_nuclei_parallel(newnuclei);
    newnuclei = new_nuclei_parallel.buildGlobalNucleiList(minDistBetweenNuclei, nuclei.size());

    // Final check to resolve overlap conflicts with existing precipitates
    std::vector<unsigned int> conflict_inds;
    safetyCheckNewNuclei(newnuclei, order_parameter_list, conflict_inds);

    newnuclei = new_nuclei_parallel.removeSubsetOfNuclei(conflict_inds, nuclei.size());

    return newnuclei;
}
