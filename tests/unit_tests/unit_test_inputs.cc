#include "../../include/userInputParameters.h"

template <int dim>
void
load_unit_test_inputs(userInputParameters<dim> &userInputs)
{
  userInputs.domain_size.push_back(1.0);
  userInputs.domain_size.push_back(1.0);
  userInputs.domain_size.push_back(1.0);
  userInputs.degree = 1;
  userInputs.subdivisions.push_back(10);
  userInputs.subdivisions.push_back(10);
  userInputs.subdivisions.push_back(10);
  userInputs.refine_factor        = 3;
  userInputs.h_adaptivity         = false;
  userInputs.max_refinement_level = 3;
  userInputs.min_refinement_level = 3;
  userInputs.skip_remeshing_steps = 1;
  userInputs.dtValue              = 1;
  userInputs.totalIncrements      = 20000;
  userInputs.finalTime            = 20.0;
  userInputs.skip_print_steps     = 1;
  userInputs.output_file_type     = "vtu";
  userInputs.nucleation_occurs    = false;
  userInputs.number_of_variables  = 1;
  userInputs.load_ICs.push_back(false);
  userInputs.load_parallel_file.push_back(false);
  userInputs.load_file_name.push_back("void");
  userInputs.load_field_name.push_back("void");
  userInputs.pp_number_of_variables = 1;
  userInputs.postProcessingRequired = false;
}
