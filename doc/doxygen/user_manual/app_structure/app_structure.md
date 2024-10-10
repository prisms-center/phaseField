# Structure of a PRISMS-PF App {#app_structure}

Each of the example application directories should contain at least eight files:

- parameters.in
- equations.cc
- ICs_and_BCs.cc
- postprocess.cc
- customPDE.h
- main.cc
- formulation.pdf
- CMakeLists.txt

It may contain two other files:
- postprocess.cc
- nucleation.cc

The following section [The Parameters File](#input_file) will discuss the structure and options for the input file, parameters.h. Users who do not need to change the governing equations or make other fundamental changes to the application will only need to modify this file.

[The App Files](#app_files) discusses the structure of the files that define an application: equations.h, ICs\_and\_BCs.h, postprocess.h, nucleation.h, customPDE.h, and main.cc. Users who want to substantially modify existing PRISMS-PF applications or create their own will find this section useful.

The final two files in each application directory are formulation.pdf and CMakeLists.txt. The file formulation.pdf describes the mathematical formulation of the problem solved by the application. Finally, the file CMakeLists.txt is used by CMake to generate the makefile. Most users will not have any reason to change CMakeLists.txt.
