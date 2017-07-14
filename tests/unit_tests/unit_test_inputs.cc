#include "../../include/userInputParameters.h"

template <int dim>
void load_unit_test_inputs(userInputParameters<dim> & userInputs){
    userInputs.domain_size.push_back(1.0);
    userInputs.domain_size.push_back(1.0);
    userInputs.domain_size.push_back(1.0);
    userInputs.degree = 1;
    userInputs.subdivisions.push_back(10);
    userInputs.subdivisions.push_back(10);
    userInputs.subdivisions.push_back(10);
    userInputs.refine_factor = 3;
    userInputs.h_adaptivity = false;
    userInputs.max_refinement_level = 3;
    userInputs.min_refinement_level = 3;
    userInputs.refine_criterion_fields.push_back(0);
    userInputs.refine_window_max.push_back(0.5);
    userInputs.refine_window_min.push_back(0.5);
    userInputs.skip_remeshing_steps = 1;
    userInputs.dtValue = 1;
    userInputs.totalIncrements = 20000;
    userInputs.finalTime = 20.0;
    userInputs.solver_type = "SolverCG";
    userInputs.abs_tol = false;
    userInputs.solver_tolerance = 1.0e-10;
    userInputs.max_solver_iterations = 1000;
    userInputs.skip_print_steps = 1;
    userInputs.output_file_type = "vtu";
    userInputs.calc_energy = false;
    userInputs.nucleation_occurs = false;
    userInputs.number_of_variables = 1;
    userInputs.load_ICs.push_back(false);
    userInputs.load_serial_file.push_back(false);
    userInputs.load_file_name.push_back("void");
    userInputs.load_field_name.push_back("void");
    userInputs.pp_number_of_variables = 1;
    userInputs.postProcessingRequired = false;

    std::vector<dealii::Tensor<2, 2*dim-1+dim/3, dealii::VectorizedArray<double> > > CIJ_list_temp;
    dealii::Tensor<2, 2*dim-1+dim/3, dealii::VectorizedArray<double> > CIJ_temp;
    CIJ_temp[0][0] = dealii::make_vectorized_array(0.0);
    CIJ_temp[1][0] = dealii::make_vectorized_array(0.0);
    CIJ_temp[2][0] = dealii::make_vectorized_array(0.0);
    CIJ_temp[0][1] = dealii::make_vectorized_array(0.0);
    CIJ_temp[1][1] = dealii::make_vectorized_array(0.0);
    CIJ_temp[2][1] = dealii::make_vectorized_array(0.0);
    CIJ_temp[0][2] = dealii::make_vectorized_array(0.0);
    CIJ_temp[1][2] = dealii::make_vectorized_array(0.0);
    CIJ_temp[2][2] = dealii::make_vectorized_array(0.0);
    CIJ_list_temp.push_back(CIJ_temp);
    userInputs.material_moduli.CIJ_list = CIJ_list_temp;


}
