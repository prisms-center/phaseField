# Version 2.0.1
Minor update to v2.0, released in November 2017. The biggest change is the introduction of a checkpoint/restart system.

Added functionality:
- A checkpoint/restart system has been added. It allows runs to be restarted if they fail or you want the run to continue for more simulated time.
- The parser to load in VTK files for initial conditions has been generalized so that it can read files generated from ParaView

Changes to the example applications:
- Fixed numerous errors in the documentation for the 'dendriticSolidification' app and edited the 'equations.h' file to use the notation from the documentation.
- Fixed typos in the documentation for the 'coupledCahnHilliardAllenCahn' app.

Bug fixes:
- Fixed the subscriptor bug that appeared at the end of simulations in debug mode.

Other changes:
- Some minor updates to the user guide, including adding a new PRISMS-PF logo.

Known issues:
- PFields only work for scalar fields and only work when the variable with index zero is a scalar field.
- Postprocessing only works for scalar fields and only when the variable with index zero is a scalar field.
- An extraneous error can appear under the following circumstances: in debug mode, with multiple cores, with adaptive meshing, using deal.II v8.4.2 or earlier. The error message includes something similar to: "Called compress(VectorOperation::insert), but the element received from a remote processor, value -2.944283206211677e-10, does not match with the value -2.944283206213795e-10 on the owner processor 60". This is a deal.II issue that they fixed in v8.5 where the tolerance for comparing numbers between processors is too tight and can be triggered by standard round-off error. It doesn't effect simulation results.

# Version 2.0:
Major update to PRISMS-PF, released in August 2017. The core library is very similar to v1.2, but the user interface is substantially changed. Most importantly, input parameters are now read from a text file, instead of via #define statements. All of the individual interface changes are not listed here, see the User Guide for the new file structure and syntax.

Added functionality:
- New option for non-uniform Dirichlet BCs (either non-uniform in space or time)
- The integral of any postprocessed field can be calculated
- Non-rectangular meshes (implemented in the new application CHiMaD_benchmark6b)

Changes to the example applications:
- New application: dendriticSolidification, solidification of a pure material, accounting for temperature evolution
- New application: anisotropyFacet, demonstrates a simple way to specify facets for strong interfacial energy anisotropy
- New applications: CHiMaD_benchmark6a and CHiMaD_benchmark6b, the new CHiMaD electrochemistry benchmark problems, with a coupled Cahn-Hilliard-Poisson system of equations

Performance improvements:
- New containers for getting variables into and residuals out of the equations.h functions yields ~15% speedup

Bug fixes:
- Unit tests now run in debug mode
- The variables that can be accessed in postprocess.h are now unrelated to those for residualRHS in equations.h

Other changes:
- Changed the syntax for accessing variables and submitting residuals in equations.h, indices are now error-checked
- To be more accurate mathematically, renamed the "ZERO_DERIVATIVE" to "NATURAL"
- User guide overhauled to use tables instead of long blocks of text, when possible
- Free energies are now treated as just another postprocessing variable
- The postprocessing.h and nucleation.h files are now optional
- Updated the README, including a link to a repository for training materials
- Version changes file changed from a general text file to markdown

Known issues:
- PFields only work for scalar fields and only work when the variable with index zero is a scalar field.
- Postprocessing only works for scalar fields and only when the variable with index zero is a scalar field.
- The formulation file for the dendriticSolidifiation application may have errors.

# Version 1.2:
Update to v1.1, released in June 2017. The core library was substantially reorganized and nucleation and postprocessing
capabilities were added.

Added functionality:
- Add postprocessing capabilities, which can be accessed in the new "postprocess.h" file in each application folder.
- New applications "nucleationModel" and "_preferential_nucleationModel" were added and include new nucleation capabilities.

Changes to the pre-built applications:
- New applications "nucleationModel" and "_preferential_nucleationModel" were added and include new nucleation capabilities.
- Changed the initial conditions for the Allen-Cahn and the Cahn-Hilliard applications to have smooth initial profiles.
- The elastic constants are now accessed by "userInputs.CIJ_List" instead of just "CIJ_List" in the functions in equations.h.
- Updated the gradient energy coefficient, the time step, and the mesh resolution of the Allen-Cahn application. Increased the
  time step for CHAC_anisotropy to its maximum stable level. Changed the mesh resolution and time step for the coupled
  Allen-Cahn/Cahn-Hilliard application. Changed the convergence tolerance, time step, and number of iterations for the precipitate
  evolution application.
- The Fickian diffusion application was overhauled to be more conceptually straightforward.
- Changed to a simpler way of specifying the initial conditions, calculating the distances by hand instead of using a pre-built,
  deal.II method.

Performance improvements:
- Removed unnecessary applications of constraints, improving run times by up to 30%.
- Updated how Dirichlet constraints are implemented, reducing the number of solver iterations for non-zero Dirichlet BCs.

Bug fixes:
- Dirichlet BCs now work for parabolic PDEs.
- Users can now load in ICs with PFields in 3D and with adaptive meshes.

Other changes:
- Restructured much of the core library. Model files were completely eliminated, pushing many of those functions into MatrixFreePDE.
  Each application now has a subclass of MatrixFreePDE called customPDE where member variables and methods can be implemented for
  specific applications.
- A wall was put in place between the define macros in the input files and the core solver code. The inputs are now housed in the
  new "userInputParameters" class.
- The core library is now linked instead of "included" in the apps. This increases the compile time the first time an app is compiled
  but decreases it substantially afterward.
- Added more in-app explanation of how to set the BCs (more comments in each ICs_and_BCs.h file)

Known issues:
- PFields only work for scalar fields and only work when the variable with index zero is a scalar field.


# Version 1.1.1:

Patch to v1.1, released in February 2017. This patch fixes some indexing bugs lurking in previous versions of the code. From a user
perspective, nothing should change (unless you had previously hit one of these bugs).

Added functionality:
- New Python script to automatically run unit tests and regression tests on the pre-built applications.
- Integration with Travis CI, for continuous integration testing.

Changes to the pre-built applications:
- Mesh adaptivity was turned on for the CHAC_anisotropy, CHAC_anisotropyRegularlized, and coupledCahnHilliardAllenCahn applications.
- Gradient coefficients were changed slightly in cahnHilliard and cahnHilliardWithAdaptivity and both applications were moved to
  second order elements with one less refinement. Mesh adaptivity parameters in cahnHilliardWithAdaptivity were tweaked.

Bug fixes:
- Fixed indexing bug that caused a seg fault if the elliptic equation variable wasn't given last.
- Fixed indexing bug that prevented the use of multiple vector variables.
- Fixed indexing bug that prevented the use of vector variables for parabolic equations.
- Fixed indexing bug that prevented the solution of multiple elliptic equations.
- Changed how input files are read, allowing PRISMS-PF to work on compilers that lack full C++11 support.

Other changes:
- Removed the model files that are no longer being used.
- Streamlined the initialization and remeshing functions.
- CMake now checks to make sure that p4est is installed.
- Fixed the description of the Cij elements for 2D anisotropic calculations.
- Updated the tensor contraction syntax for deal.II 8.4.

Known issues:
- None


# Version 1.1:

Update to v1.0, released in January 2017.

Added functionality:
- Periodic boundary conditions are now available. Choose “PERIODIC” as the BC type in the function setBCs() in ICs_and_BCs.h.
- PRISMS-PF now automatically constrains the solution if only periodic/Neumann BCs are used for a component of a variable by pinning
  the value of that variable to zero at the origin.
- Added more options for controlling the spacing between outputs. Users can now choose between equal spacing, log spacing, and a specified
  number of outputs per decade or they can input a list of time steps at which to output.
- Added the ability to control how often the progress of a run prints to the screen. Set “skipPrintSteps” in parameters.h.
- Ability to change the output file format between “.vtu” and “.vtk”. Used in the application allenCahn_pfield.
- When mesh adaptivity is turned on, the mesh will adapt before the first time step.

Changes to the pre-built applications:
- New application: singlePrecipitateKKS is a 3D simulation of 1/8 of a precipitate. Governing equations are similar to precipitateEvolution
  but use the KKS model instead of the WBM model.
- New application: allenCahn_pfield uses the PField capability of PRISMS IntegrationTools to load the initial condition from .vtk files
- New application: precipitateEvolution_pfunction uses the PFunction capability of PRISMS IntegrationTools to load parameters (both
  and functions) into the code.
- New applications: CHiMaD_benchmark1a and CHiMaD_benchmark2a are phase field benchmark problems for spinodal decomposition and Ostwald
  ripening, respectively. More about the benchmark problems can be found at: https://pages.nist.gov/chimad-phase-field/
- Added formulation files for the CHAC_anisotropy and CHAC_anisotropyRegularized applications.
- Boundary conditions for grainGrowth changed to periodic.
- Boundary condition for the mechanics application changed for higher deformation.
- Moduli changed in precipitateEvolution so that the evolution is dominated by the misfit, not the interfacial energy.

Bug fixes:
- Fixed constraint bug for purely elliptic problems. The solution in the domain was correct, but the boundary value for Dirichlet BCs was
  double what it was supposed to be.
- Removed unnecessary LAPACK check in the CMakeLists.txt files.

Other changes:
- Installation instructions updated, both in the README and in the user guide.
- The “init” function was split into an “init” function and a “reinit” function. The latter is used for reinitializing the system for
  adaptive meshing.
- Removed call of “norm_sqr”, which no longer exists in deal.II v8.5.

Known issues:
- For some compilers that lack full C++11 support (including the Intel compilers on the Stampede cluster), errors are generated because
  we initialize some vectors with C-style arrays. This will be fixed in the near future.

# Version 1.0:

This is the first release version of PRISMS-PF, released in August 2016. The code has been substantially reworked since version 0.9.3.

Major changes include:
- New user interface
- New user guide
- Mesh adaptivity capability
- Many new applications
- Improved tests, including tests versus analytical results


# Version 0.9.3:

Patch to version 0.9.2, released in July 2016. Most changes are specific to coupled Cahn-Hilliard-Allen-Cahn-Mechanics calculations.

Bug fixes:
- Made necessary fixes for compatibility Deal.II v8.4.

Added functionality:
- Support for concentration-dependent stress-free transformation strains for coupled Cahn-Hilliard-Allen-Cahn-Mechanics calculations.
- Support for heterogeneous mechanics for coupled Cahn-Hilliard-Allen-Cahn-Mechanics calculations.
- Option to shift the concentration by a constant value to reach a desired average concentration for coupled
  Cahn-Hilliard-Allen-Cahn-Mechanics calculations.
- Free energy integrator for coupled Cahn-Hilliard-Allen-Cahn-Mechanics calculations.
- Option to use an absolute tolerance for mechanics calculations.
- New speed comparison test versus a finite difference code for coupled Cahn-Hilliard-Allen-Cahn calculations.
- New application "bPPE_pfunction" that is a mirror of "betaPrimePrecipitateEvolution", but shows how PFunctions from the PRISMS
  IntegrationTools can be used to input parameters.

Performance improvements:
- Added option for fewer than 3 structural order parameters for coupled Cahn-Hilliard-Allen-Cahn-Mechanics calculations, which can
  cut run time in half.
- Refactored mechanics functions for coupled Cahn-Hilliard-Allen-Cahn-Mechanics calculations, making use of built-in Deal.II tensor
  functions. This improves speed and reduces code.

Other changes:
- Restructured input files to separate numerical parameters, ICs/BCs, and the residual equations. This is a stepping stone to a more
  substantial interface overhaul in the next release.
- Updated/tweaked the example applications to make sure interfaces were well resolved, etc.

Known issues:
- The PFunctions in the application "bPPE_pfunction" can currently only be used for constants. The PFunctions only return doubles and
  thus cannot be used for functions like the free energy that must be vectorized arrays. This issue is being actively worked on.


# Version 0.9.2:

Patch to version 0.9.1, released in January 2016.

Bug fixes:
- In versions 0.9 and 0.9.1, unlike version 0.8, the fields in the calculation were updated sequentially, with the updated value being
  used for other calculations in that time step (e.g. the concentration is updated using the Cahn-Hilliard equation, then that updated
  value is used when updating the order parameter via the Allen-Cahn equation). In version 0.8, all of the updates to the fields used
  the value of the fields from the previous time step. This was the source of the different solutions between the versions. We
  reverted to the approach taken in version 0.8.

Added functionality:
- A new test suite has been added. The tests can be found in the "tests" directory. Tests include a regression test, testing that the
  solutions exhibit the expected order of accuracy in both time and space, a comparison against a finite difference code, and tests
  comparing the accuracy of first, second, and third order elements at a range of element sizes.

Performance improvements:
- Refactoring of "computeRHS" has led to a substantial performance increase for solving parabolic equations. Most of the improved
  performance came from no longer storing field metadata in std::map containers.

Known issues:
- None


# Version 0.9.1:

Patch to version 0.9, released in December 2015. This patch fixes a number of bugs that were introduced between versions 0.8 and 0.9.

Bug fixes:
- Fixed error in strain calculation in getRHS, previously it missed some terms in 3D calculations
- Fixed formatting for output filenames after the 999,999th output
- Fixed a sign error in the residuals in the beta prime precipitate application
- Fixed symmetry error in generation of the stiffness tensor for 2D anisotropic calculations
- Fixed error in the mechanics solver where solution was erroneously set to the change in the solution
- Updated the namespace for MPI calls to comply with Deal.II v8.3 specifications

Added functionality:
- Changed mesh generation to allow for elements with a different aspect ratio than the domain

Known issues:
- Code runs more slowly than version 0.8, possibly due to compiler flags set by Deal.II that are disabling vectorization
- Yields different solutions for time evolution equations than version 0.8


# Version 0.9:

Released in September 2015 before the PRISMS Workshop. A large portion of the code was refactored so that much more of the source code is
shared between applications. In the previous version, each application was nearly independent.

Known issues:
- Code runs more slowly than version 0.8, possibly due to compiler flags set by Deal.II that are disabling vectorization
- Yields different solutions for time evolution equations than version 0.8


# Version 0.8:

First public release of the code.
