PRISMS-PF
=================
<B>Code repository:</B> https://github.com/prisms-center/phaseField <br>
<B>Code documentation:</B> https://goo.gl/00y23N <br>
<B>User registration link:</B> http://goo.gl/forms/GXo7Im8p2Y

<B>Version information:</B>

This version of the code, 1.0, is the first release version of PRISMS-PF. This version provides many updates to the pre-release versions. For more information concerning the differences between versions, please consult version_changes.txt.   

<B>What is PRISMS-PF?</B>

  PRISMS-PF is a high-performance Finite Element Method (FEM) code for phase field modeling simulations and solving related systems of partial differential equations. It has a flexible user interface designed to allow users to easily modify existing applications as well as developing their own. In terms of performance, PRISMS-PF has been shown to be competitive with state-of-the-art finite difference codes. PRISMS-PF has been tested on thousands of processors and has been used in simulations with over one billion degrees of freedom.
  
  This code is developed by the PRedictive Integrated Structural
  Materials Science (PRISMS) Center [http://www.prisms-center.org/]
  at University of Michigan which is supported by the U.S. Department 
  of Energy (DOE), Office of Basic Energy Sciences, Division of Materials Sciences 
  and Engineering under Award #DE-SC0008637

<B>Quick Start Guide:</B>

For detailed instructions on how to download and use PRISMS-PF, please consult the PRISMS-PF Users Guide (the file prismspf_users_guide.pdf). An abbreviated version of the instructions is given below.

<I>Installation:</I> 

1) Install deal.II (version 8.4.1 recommended)<br>
  + Download CMake [http://www.cmake.org/download/]
  + Add CMake to your path (e.g. $ PATH="/path/to/cmake/Contents/bin":"$PATH"), preferably in a shell configuration file 
  + Download Deal.II binaries (OSX and Linux) or  Virtual Machine (VMI) from https://www.dealii.org/download.html 
  + If a Deal.II binary is downloaded, open it and follow the instructions in the terminal window <br>

2) Clone the PRISMS-PF GitHub repo https://github.com/prisms-center/phaseField<br>
  + $ git clone https://github.com/prisms-center/phaseField.git <br>
  + $ cd phaseField <br>
  + $ git checkout master <br>
  
<I>Updates:</I> 

Since PRISMS-PF is still under active development,
  regular code and documentation updates are pushed to the upstream
  repo (https://github.com/prisms-center/phaseField) and we strongly
  recommend users to synchronize their respective clones/forks at regular
  intervals or when requested by the developers through the
  announcements on the mailing list. 

<I>Running a Pre-Built Application:</I> 

  Entering the following commands will run one of the pre-built example applications (the Cahn-Hilliard spinodal decomposition application in this case):<br> 
  + $ cd applications/cahnHilliard <br>
  For debug mode [default mode, very slow]: <br>
  + $ cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Debug <br>
  For optimized mode:<br>
  + $ cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Release <br>
  and <br>
  + $ make <br><br>
  Execution (serial runs): <br>
  + $ make run <br>
  Execution (parallel runs): <br>
  + $ mpirun -np nprocs ./main <br>
  [here nprocs denotes the number of processors]
  
<I>Visualization:</I> 

  Output of the primal fields and postprocessed fields is in standard vtk 
  format (parallel:*.pvtu, serial:*.vtu files) which can be visualized with the 
  following open source applications:
  1. VisIt (https://wci.llnl.gov/simulation/computer-codes/visit/downloads)
  2. Paraview (http://www.paraview.org/download/)

<I>Getting started:</I>

  Examples of various initial boundary value problems (IBVP's) are located under the 
  applications/ folder. The easiest way to get started on the code is to 
  run the example applications in this folder.

  THe example applications are intended to serve as (1) Demonstration of the
  capabilities of this library, (2) Provide a framework for
  further development of specialized/advanced applications by
  users. 

  Applications that are still under development/testing are preceded by an
  underscore. 

<B>Documentation:</B>

  The PRISMS-PF Users Guide provides extensive documentation on installing PRISMS-PF, running and visualizing simulations, and the structure of the input files.
  
  Doxygen-generated documentation can be viewed in one of two ways: 
  + Open html/index.html in any web browser <br>
  (OR)<br>
  + https://htmlpreview.github.io/?https://raw.githubusercontent.com/prisms-center/phaseField/master/html/index.html
 	
<B>License:</B>

  GNU Lesser General Public License (LGPL). Please see the file
  LICENSE for details.

<B>Mailing Lists:</B>
  
 + prismsphaseField.users@umich.edu	
 + prismsphaseField.dev@umich.edu  

<B>Further information, questions, issues and bugs:</B>

  Contact the developers at prismsphaseField.dev@umich.edu  



