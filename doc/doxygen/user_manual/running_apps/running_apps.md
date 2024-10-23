# Running a PRISMS-PF Example App {#running_apps}
## The Example Apps

After deal.II and PRISMS-PF are downloaded, you can run the pre-built PRISMS-PF example applications. At this time, the example applications include:
- allenCahn: An implementation of the Allen-Cahn equation for two phases. (2D)
- cahnHilliard: An implementation of the Cahn-Hilliard equation for two phases. (2D)
- coupledCahnHilliardAllenCahn: An implementation of the coupled Cahn-Hilliard/Allen-Cahn set of equations. (2D)
- CHAC_performance_test: An implementation of the coupled Cahn-Hilliard/Allen-Cahn set of equations for two growing particles used for benchmarking purposes. (3D)
- CHAC\_anisotropy: Coupled Cahn-Hilliard/Allen-Cahn equations with weakly anisotropic interfacial energy. (2D)
- CHAC\_anisotropyRegularized: Like CHAC\_anisotropy, but with a regularization term to permit strongly anisotropic interfacial energy. (2D)
- anisotropyFacet: A different strong anisotropy formation than in CHAC\_anisotropyRegularized that is easier to specify particular facets in the Wulff shape. (2D)
- steadyStateAllenCahn: An implementation of coupled Allen-Cahn and steady-state Allen-Cahn equations as a demonstration of the nonlinear solver. (2D)
- fickianDiffusion: An implementation of the diffusion equation with a time-dependent source term. (2D)
- mechanics: An implementation of linear elasticity for a material in uniaxial tension. (3D)
- CHiMaD\_benchmark1a: An implementation of the CHiMaD spinodal decomposition benchmark problem. (2D)
- CHiMaD\_benchmark2a: An implementation of the CHiMaD Ostwald ripening benchmark problem. (2D)
- CHiMaD\_benchmark3: An implementation of the CHiMaD dendritic solidification benchmark problem. (2D)
- CHiMaD\_benchmark6a: An implementation of the CHiMaD electrochemistry benchmark problem. (2D)
- CHiMaD\_benchmark6b: An implementation of the CHiMaD electrochemistry benchmark problem with a curved domain. (2D)
- dendriticSolidification: An implementation of a solidification model for a pure material resulting in the growth of a dendrite. (2D)
- eshelbyInclusion: An implementation of linear elasticity for a spherical inclusion. (3D)
- grainGrowth: An implementation of coupled Allen-Cahn equations simulating grain growth in two dimensions. (2D)
- grainGrowth_dream3D: An implementation of coupled Allen-Cahn equations simulating grain growth in two dimensions with an initial microstructure imported from Dream3D. (2D)
- precipiateEvolution: An implementation of the coupled Cahn-Hilliard/Allen-Cahn/Linear Elasticity equations often used in phase field simulation of precipitate evolution. (2D)
- precipiateEvolution_pfunction: Like precipitateEvolution, but loads inputs using PRISMS IntegrationTools. (2D)
- MgNd_precipitate_single_Bppp: Similar to precipiateEvolution, but uses the KKS model rather than the WBM model for the free energy functional. The parameters are set for \f$\beta\f$''' precipitates in an Mg-Nd alloy from [this publication](https://www.sciencedirect.com/science/article/pii/S1359645417305281). (3D)
- nucleationModel: KKS precipitation model that makes use of the PRISMS-PF explicit nucleation capabilities. (2D)
- nucleationModel_preferential: Like nucleationModel, but with a zone with an increased nucleation rate to simulate a grain boundary. (2D)

A directory for each of these apps can be found in the [applications directory](https://github.com/prisms-center/phaseField/tree/master/applications) (i.e. phaseField/applications). The apps contain a formulation file giving the governing equations. In addition to the 24 apps listed above, some app names may be preceded by an underscore. The underscore is used to denote apps that are still under active development.

## Compiling and Running the Allen-Cahn Example App

From the ''phaseField'' directory one can run the Allen-Cahn example application through to following terminal commands:
```
$ cd applications/allenCahn/
$ cmake .
$ make debug
$ mpirun -n 1 main
```
The first command moves from the ''phaseField'' directory to the directory of the Allen-Cahn example. The second command checks that core PRISMS-PF library has been compiled, (re-)compiles it if necessary, and creates a \emph{makefile} using CMake. The third command compiles the executable in ''debug'' mode, which enables a number of exception checks in the code and adds debugging information that can be used by a debugger (e.g. gdb). The fourth command runs the program using a single processor.

As the program runs, information from each time step outputs to the terminal window. After the simulation is complete, a summary the time taken in a few major sections of the code and the total wall time is printed to the terminal window.

Here is a screenshot of typical output from CMake as you create the \emph{makefile}:
![](cmake_output_v2.png)

Don't worry if the output isn't exactly the same as what you see, the details of some of the messages depend on your operating system and which compilers you have installed. The important part is that the bottom three messages are ''Configuring done'', ''Generating done'', and ''Build files have been written to: ...''. In the future, entering ''\$ cmake .'' will result in a shorter set of messages because CMake caches some variables from the last time it was run. As a result, you can omit the CMake step for future simulations as long as the path name to your current directory is unchanged and your installation of deal.II is unchanged.

Here is a screenshot of typical output from the compiler as you compile the executable:
![](compile_output.png)

Depending on your version of deal.II, different warnings may appear as you compile. Common warnings include the use of functions that deal.II has marked as depricated (as in the screenshot above) and unused type definitions. In this case, PRISMS-PF uses these functions for backward compatability with deal.II version 8.4.x. We will switch to the updated functions in the near future.

Once the simultation is complete, the terminal output at the end of the simulation should look like:

![](allenCahn_output.png)

## What Can Go Wrong
If you were able to enter all of the commands in the previous section and get output similar to the screenshots, congratulations! you just ran your first PRISMS-PF simulation. If not, you may be experiencing one of the common issues listed below.

If CMake gives an error message like this:
![](cmake_no_period.png)

Then you likely forgot the period at the end of the command ```$ cmake .```.

If CMake gives an error message like this:
![](cmake_no_dealii.png)

CMake cannot find your installation of deal.ii. This issue is probably caused by the lack of an environment variable pointing to the directory containing your deal.II library. You can check this with the following command:
```
$ echo $DEAL_II_DIR
```
The terminal should then output the path to the deal.II library. For example in Mac OS, the deal.II directory may be ''/Applications/deal.II.app/Resources''. If  DEAL_II_DIR contains a path, go there to see if deal.II is actually installed there. If DEAL\_II\_DIR is incorrect or empty, you should set it to the directory of the correct path of your deal.II installation with the following command:
```
$ export DEAL_II_DIR=/path/to/dealii
```
Environment variables are erased when you close your shell. To have this path set every time you open a shell (i.e. every time you open a new terminal window), you can add the command above to your shell profile (e.g. .bashrc if you use a bash shell). If you are still having problems, there may be an issue with your deal.II installation. Please consult the deal.II website for instructions.

If CMake cannot successfully detect a C++ compiler, it will generate an error message. The most common cause for this is that the machine runs the Mac OS operating system with an outdated version of the Clang compiler. Upgrading your OS to version 10.11 or newer, updating Xcode, and (re)installing the Xcode command line tools may help. Alternatively, you can install a certain version of the deal.II package that was developed to sidestep this issue:
https://github.com/dealii/dealii/releases/download/v8.3.0/dealii-8.3.0.nocxx14.dmg

Most of issues users have had are during the CMake step. If the fixes suggested above don't work for your or you have an issue not covered by this list, please contact the PRISMS-PF users list: prisms-pf-users@googlegroups.com. If you are not already on the list, please submit a join request through Google Groups or send an email with ''SUBSCRIBE'' in the subject line to prismsphasefield.dev@umich.edu. As users come across new issues, we will add them (and suggested fixes) to this section.

## Visualizing the Results of the Simulation
Once you have successfully run a simulation, you will likely want to visualize the results.  PRISMS-PF output files are generated in the popular VTK format, as a series of .vtu and .pvtu files. Two common open-soure, multi-platform visualization tools for these types of files are VisIt and ParaView. Instructions for downloading this software can be found at their respective websites:

VisIt: https://wci.llnl.gov/simulation/computer-codes/visit/
\
ParaView: http://www.paraview.org/

To get you started, here is a brief tutorial on how to use VisIt to visualize your simulation results. For more detailed instructions, please consult the [VisIt manual](https://wci.llnl.gov/simulation/computer-codes/visit/manuals).

After launching the VisIt application, click ''Open'' and find the directory for the Allen Cahn example:
![](visit_open.png)

and select ''solution-\*.pvtu'':
![](visit_path.png)

Next, click ''Add'', hover over ''Pseudocolor'', and select ''n'':
![](visit_add.png)

Next, click ''Draw'' to make a plot:
![](visit_pseudocolor.png)

The VisIt window will now have a plot of the initial condition of the order parameter:
![](visit_result.png)

The result at other time steps can be visualized by dragging the time bar, or using the controls directly below the time bar. Dragging the time bar to the end will display the final result of the simulation:
![](visit_time_bar.png)
![](visit_final_result.png)

VisIt has a wide variety of capabilities for visualizing 1D, 2D, and 3D data, including postprocessing of output fields (e.g. to obtain the result of mathematical expressions involving one or more output field). VisIt also has a powerful Python interface to provide scripting capabilities. More more instructions on how to use these, and other, features of VisIt, please consult the [VisIt manual](https://wci.llnl.gov/simulation/computer-codes/visit/manuals).

## Running the Other Example Applications
Running the other example applications is as simple as going to the directory for the application of interest and repeating the steps in [the Allen-Cahn section above](#Compiling-and-Running-the-Allen-Cahn-Example-App). For example, to run the Cahn-Hilliard example application, when you are currently in the directory for the Allen-Cahn example application, you would type the following commands:
```
$ cd ../cahnHilliard/
$ cmake .
$ make debug
$ mpirun -n 1 main
```

## Running in Release Mode and with Multple Processors
Once you are sure that your code works as expected, you can compile it in ''Release Mode'', which is approximately 10x faster than ''Debug Mode''. However, no exceptions are checked in Release Mode and therefore the code may compile and run even if it contains several errors. Therefore, we strongly recommend that all new code is tested in Debug Mode before switching to Release Mode. To compile in Release mode, replace ```$ make debug``` with ```$ make release```.

One of the strengths of PRISMS-PF is that it can be run in parallel with almost no extra effort from the user. To run a simulation in parallel, replace the ''1'' in the ```mpirun``` command with the desired number of processor cores. The deal.II packages on their website already contain the MPI library, so no extra software has to be downloaded to use multiple cores.

From the directory of an example application, a simulation can be run in Release Mode on 4 cores using the following commands:
```
$ cmake .
$ make release
$ mpirun -n 4 main
```
