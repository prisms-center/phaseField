// Method for the userInputParameters class
#include "../include/inputFileReader.h"

inputFileReader::inputFileReader(dealii::ParameterHandler & parameter_handler, std::string input_file_name, unsigned int number_of_variables, unsigned int number_of_materials, unsigned int number_of_pp_variables){

    // Declare all of the entries
    parameter_handler.declare_entry("Number of dimensions","-1",dealii::Patterns::Integer(),"The number of dimensions for the simulation.");

    parameter_handler.declare_entry("Domain size X","-1",dealii::Patterns::Double(),"The size of the domain in the x direction.");
    parameter_handler.declare_entry("Domain size Y","-1",dealii::Patterns::Double(),"The size of the domain in the y direction.");
    parameter_handler.declare_entry("Domain size Z","-1",dealii::Patterns::Double(),"The size of the domain in the z direction.");
    parameter_handler.declare_entry("Element degree","1",dealii::Patterns::Integer(),"The polynomial order of the finte element.");
    parameter_handler.declare_entry("Subdivisions X","1",dealii::Patterns::Integer(),"The number of mesh subdivisions in the x direction.");
    parameter_handler.declare_entry("Subdivisions Y","1",dealii::Patterns::Integer(),"The number of mesh subdivisions in the y direction.");
    parameter_handler.declare_entry("Subdivisions Z","1",dealii::Patterns::Integer(),"The number of mesh subdivisions in the z direction.");
    parameter_handler.declare_entry("Refine factor","-1",dealii::Patterns::Integer(),"The number of initial refinements of the coarse mesh.");

    parameter_handler.declare_entry("Mesh adaptivity","false",dealii::Patterns::Bool(),"Whether to enable mesh adaptivity.");
    parameter_handler.declare_entry("Max refinement level","-1",dealii::Patterns::Integer(),"The maximum level of refinement.");
    parameter_handler.declare_entry("Min refinement level","-1",dealii::Patterns::Integer(),"The minimum level of refinement.");
    parameter_handler.declare_entry("Refinement criteria fields","0",dealii::Patterns::List(dealii::Patterns::Integer()),"The list of fields used to determine mesh refinement.");
    parameter_handler.declare_entry("Refinement window max","",dealii::Patterns::List(dealii::Patterns::Anything()),"The upper limit for refinement for each of the criteria fields.");
    parameter_handler.declare_entry("Refinement window min","",dealii::Patterns::List(dealii::Patterns::Anything()),"The lower limit for refinement for each of the criteria fields.");
    parameter_handler.declare_entry("Steps between remeshing operations","1",dealii::Patterns::Integer(),"The number of time steps between mesh refinement operations.");

    parameter_handler.declare_entry("Number of time steps","-1",dealii::Patterns::Integer(),"The time step size for the simulation.");
    parameter_handler.declare_entry("Time step","-0.1",dealii::Patterns::Double(),"The time step size for the simulation.");
    parameter_handler.declare_entry("Simulation end time","-0.1",dealii::Patterns::Double(),"The value of simulated time where the simulation ends.");

    parameter_handler.declare_entry("Linear solver","SolverCG",dealii::Patterns::Anything(),"The linear solver (currently only SolverCG).");
    parameter_handler.declare_entry("Use absolute convergence tolerance","false",dealii::Patterns::Bool(),"Whether to use an absolute tolerance for the linear solver (versus a relative tolerance).");
    parameter_handler.declare_entry("Solver tolerance value","1.0e-3",dealii::Patterns::Double(),"The tolerance for the linear solver (either absolute or relative).");
    parameter_handler.declare_entry("Maximum allowed solver iterations","10000",dealii::Patterns::Integer(),"The maximum allowed number of iterations the linear solver is given to converge before being forced to exit.");

    parameter_handler.declare_entry("Write output","true",dealii::Patterns::Bool(),"Whether or not to write output files for the simulation.");
    parameter_handler.declare_entry("Output file type","vtu",dealii::Patterns::Anything(),"The output file type (either vtu or vtk).");
    parameter_handler.declare_entry("Output condition","EQUAL_SPACING",dealii::Patterns::Anything(),"The time step size for the simulation.");
    parameter_handler.declare_entry("List of time steps to output","0",dealii::Patterns::Anything(),"The list of time steps to output, used for the LIST type.");
    parameter_handler.declare_entry("Number of outputs","10",dealii::Patterns::Integer(),"The number of outputs (or number of outputs per decade for the N_PER_DECADE type).");
    parameter_handler.declare_entry("Skip print steps","1",dealii::Patterns::Integer(),"The number of time steps between updates to the screen.");
    parameter_handler.declare_entry("Calculate the free energy","false",dealii::Patterns::Bool(),"Whether to calculate the integrated free energy.");

    parameter_handler.declare_entry("Allow nucleation","false",dealii::Patterns::Bool(),"Whether to enable the explicit nucleation capabilties.");

    // Declare entries regarding the governing equations
    parameter_handler.declare_entry("Number of equations","-1",dealii::Patterns::Integer(),"The number of governing equations being solved.");

    for (unsigned int i=0; i<number_of_variables; i++){
        std::string equation_text = "Equation ";
        equation_text.append(dealii::Utilities::int_to_string(i));

        parameter_handler.enter_subsection(equation_text);
        {
            parameter_handler.declare_entry("Variable name","var",dealii::Patterns::Anything(),"The name of the field variable.");
            parameter_handler.declare_entry("Variable type","SCALAR",dealii::Patterns::Anything(),"Whether the variable is a SCALAR or a VECTOR.");
            parameter_handler.declare_entry("Equation type","PARABOLIC",dealii::Patterns::Anything(),"Whether the governing equation is PARABOLIC or ELLIPTIC.");

            parameter_handler.declare_entry("Need variable value","true",dealii::Patterns::Bool(),"Whether the value of the variable is needed for any of the residual equations.");
            parameter_handler.declare_entry("Need variable gradient","true",dealii::Patterns::Bool(),"Whether the gradient of the variable is needed for any of the residual equations.");
            parameter_handler.declare_entry("Need variable hessian","false",dealii::Patterns::Bool(),"Whether the hessian of the variable is needed for any of the residual equations.");

            parameter_handler.declare_entry("Need value residual term","true",dealii::Patterns::Bool(),"Whether the residual equation has a term proportional to the value of the test function.");
            parameter_handler.declare_entry("Need gradient residual term","true",dealii::Patterns::Bool(),"Whether the residual equation has a term proportional to the gradient of the test function.");

            parameter_handler.declare_entry("Need variable value (LHS)","false",dealii::Patterns::Bool(),"Whether the value of the variable is needed for any of the LHS residual equations.");
            parameter_handler.declare_entry("Need variable gradient (LHS)","false",dealii::Patterns::Bool(),"Whether the gradient of the variable is needed for any of the LHS residual equations.");
            parameter_handler.declare_entry("Need variable hessian (LHS)","false",dealii::Patterns::Bool(),"Whether the hessian of the variable is needed for any of the LHS residual equations.");

            parameter_handler.declare_entry("Need value residual term (LHS)","false",dealii::Patterns::Bool(),"Whether the LHS residual equation has a term proportional to the value of the test function.");
            parameter_handler.declare_entry("Need gradient residual term (LHS)","false",dealii::Patterns::Bool(),"Whether the LHS residual equation has a term proportional to the gradient of the test function.");
        }
        parameter_handler.leave_subsection();
    }

    // Declare entries regarding the elastic constants
    parameter_handler.declare_entry("Number of materials","0",dealii::Patterns::Integer(),"The number of materials with differing elastic constants.");


    for (unsigned int i=0; i<number_of_materials; i++){
        std::string material_text = "Material ";
        material_text.append(dealii::Utilities::int_to_string(i));

        parameter_handler.enter_subsection(material_text);
        {
            parameter_handler.declare_entry("Elastic constants","-1.0",dealii::Patterns::List(dealii::Patterns::Double()),"The elastic constants for the material (see the user guide for the ordering).");
            parameter_handler.declare_entry("Material symmetry","ISOTROPIC",dealii::Patterns::Anything(),"The elastic symmetry of the material (see the user guide for the list of options).");
        }
        parameter_handler.leave_subsection();
    }

    // Declare entries for reading initial conditions from file
    parameter_handler.declare_entry("Load initial conditions","void",dealii::Patterns::Anything(),"Whether to load the initial conditions for each variable from file.");
    parameter_handler.declare_entry("Load serial file","void",dealii::Patterns::Anything(),"Whether all processors should read from a single file (versus each reading from separate files).");
    parameter_handler.declare_entry("File names","void",dealii::Patterns::Anything(),"The file name to load from for each variable.");
    parameter_handler.declare_entry("Variable names in the files","void",dealii::Patterns::Anything(),"What each variable is named in the field being loaded.");

    // Declare entries for the postprocessing variables
    parameter_handler.declare_entry("Number of postprocessing variables","-1",dealii::Patterns::Integer(),"The number of governing equations being solved.");

    for (unsigned int i=0; i<number_of_pp_variables; i++){
        std::string pp_var_text = "Postprocessing variable ";
        pp_var_text.append(dealii::Utilities::int_to_string(i));

        parameter_handler.enter_subsection(pp_var_text);
        {
            parameter_handler.declare_entry("Variable name","var",dealii::Patterns::Anything(),"The name of the field variable.");
            parameter_handler.declare_entry("Variable type","SCALAR",dealii::Patterns::Anything(),"Whether the variable is a SCALAR or a VECTOR.");

            parameter_handler.declare_entry("Need variable value","true",dealii::Patterns::Bool(),"Whether the value of the variable is needed for any of the residual equations.");
            parameter_handler.declare_entry("Need variable gradient","true",dealii::Patterns::Bool(),"Whether the gradient of the variable is needed for any of the residual equations.");
            parameter_handler.declare_entry("Need variable hessian","false",dealii::Patterns::Bool(),"Whether the hessian of the variable is needed for any of the residual equations.");

            parameter_handler.declare_entry("Need value residual term","true",dealii::Patterns::Bool(),"Whether the residual equation has a term proportional to the value of the test function.");
            parameter_handler.declare_entry("Need gradient residual term","true",dealii::Patterns::Bool(),"Whether the residual equation has a term proportional to the gradient of the test function.");
        }
        parameter_handler.leave_subsection();
    }


    // Read from the input file
    parameter_handler.read_input(input_file_name);

}
