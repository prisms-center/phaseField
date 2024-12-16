var manual =
[
    [ "PRISMS-PF Overview", "manual.html#autotoc_md15", null ],
    [ "Structure of the PRISMS-PF User Guide", "manual.html#autotoc_md16", null ],
    [ "Installation", "install.html", "install" ],
    [ "Running a PRISMS-PF Example App", "running_apps.html", [
      [ "The Example Apps", "running_apps.html#autotoc_md55", null ],
      [ "Compiling and Running the Allen-Cahn Example App", "running_apps.html#autotoc_md56", null ],
      [ "What Can Go Wrong", "running_apps.html#autotoc_md57", null ],
      [ "Visualizing the Results of the Simulation", "running_apps.html#autotoc_md58", null ],
      [ "Running the Other Example Applications", "running_apps.html#autotoc_md59", null ],
      [ "Running in Release Mode and with Multple Processors", "running_apps.html#autotoc_md60", null ]
    ] ],
    [ "Structure of a PRISMS-PF App", "app_structure.html", null ],
    [ "The Input File: parameters.in", "input_file.html", [
      [ "General Structure of the Input File", "input_file.html#autotoc_md17", null ],
      [ "Tables Describing Each Block of Input Parameters", "input_file.html#autotoc_md18", [
        [ "Dimensionality", "input_file.html#autotoc_md19", null ],
        [ "Computational Domain", "input_file.html#autotoc_md20", null ],
        [ "Element Parameters", "input_file.html#autotoc_md21", null ],
        [ "Mesh Adaptivity (optional)", "input_file.html#autotoc_md22", null ],
        [ "Time Stepping", "input_file.html#autotoc_md23", null ],
        [ "Linear Solver Parameters for Each TIME_INDEPENDENT Equation (optional, see Note 1 below for details)", "input_file.html#autotoc_md24", null ],
        [ "Shared Nonlinear Solver Parameters (optional, see Note 2 below for details)", "input_file.html#autotoc_md25", null ],
        [ "Nonlinear Solver Parameters for Each Nonlinear Equation (optional, see Note 2 below for details)", "input_file.html#autotoc_md26", null ],
        [ "Output (optional)", "input_file.html#autotoc_md27", null ],
        [ "Checkpoints (optional, see Note 3 below for details)", "input_file.html#autotoc_md28", null ],
        [ "Boundary Conditions (see Note 4 below for details)", "input_file.html#autotoc_md29", null ],
        [ "Loading Initial Conditions from File (optional, see Note 5 below for details)", "input_file.html#autotoc_md30", null ],
        [ "Shared Nucleation Parameters (optional, see Note 6 below for details)", "input_file.html#autotoc_md31", null ],
        [ "Nucleation Parameters for Each Nucleating Variable (optional, see Note 6 below for details)", "input_file.html#autotoc_md32", null ],
        [ "Grain Remapping Parameters (optional, see Note 7 below for details)", "input_file.html#autotoc_md33", null ],
        [ "Model constants (optional, see Note 8 below for details)", "input_file.html#autotoc_md34", null ],
        [ "Note 1: Linear Solver Parameters", "input_file.html#autotoc_md35", null ],
        [ "Note 2: Nonlinear Solver Parameters", "input_file.html#autotoc_md36", null ],
        [ "Note 3: Checkpoint/Restart", "input_file.html#autotoc_md37", null ],
        [ "Note 4: Boundary Conditions", "input_file.html#autotoc_md38", null ],
        [ "Note 5: Loading Initial Conditions from File", "input_file.html#autotoc_md39", null ],
        [ "Note 6: Nucleation", "input_file.html#autotoc_md40", null ],
        [ "Note 7: Grain Remapping Parameters", "input_file.html#autotoc_md41", null ],
        [ "Note 8: Model Constants", "input_file.html#autotoc_md42", null ]
      ] ]
    ] ],
    [ "Other App Files", "app_files.html", [
      [ "equations.cc", "app_files.html#autotoc_md0", [
        [ "loadVariableAttributes", "app_files.html#autotoc_md1", null ],
        [ "explicitEquationRHS", "app_files.html#autotoc_md2", [
          [ "A note on types", "app_files.html#autotoc_md3", null ]
        ] ],
        [ "nonExplicitEquationRHS", "app_files.html#autotoc_md4", null ],
        [ "equationLHS", "app_files.html#autotoc_md5", null ]
      ] ],
      [ "ICs_and_BCs.cc", "app_files.html#autotoc_md6", [
        [ "setInitialCondition", "app_files.html#autotoc_md7", null ],
        [ "setNonUniformDirichletBCs", "app_files.html#autotoc_md8", null ]
      ] ],
      [ "postprocess.cc", "app_files.html#autotoc_md9", [
        [ "loadPostProcessorVariableAttributes", "app_files.html#autotoc_md10", null ],
        [ "postProcessedFields", "app_files.html#autotoc_md11", null ]
      ] ],
      [ "nucleation.cc", "app_files.html#autotoc_md12", null ],
      [ "customPDE.h", "app_files.html#autotoc_md13", null ],
      [ "main.cc", "app_files.html#autotoc_md14", null ]
    ] ],
    [ "Tests", "tests.html", null ],
    [ "Creating Custom PRISMS-PF Apps", "custom_apps.html", null ]
];