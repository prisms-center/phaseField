// Method for the userInputParameters class
#include "../include/inputFileReader.h"

inputFileReader::inputFileReader(dealii::ParameterHandler & parameter_handler, std::string input_file_name){

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

    // Read from the input file
    parameter_handler.read_input(input_file_name);

}
