// Parameter list for the precipitate evolution example application
// All strictly numerical parameters should be set in this file

// =================================================================================
// Set the number of dimensions (1, 2, or 3 for a 1D, 2D, or 3D calculation)
// =================================================================================
#define problemDIM 2

// =================================================================================
// Set the length of the domain in all three dimensions
// =================================================================================
// Each axes spans from zero to the specified length
#define scaleFactor 2.0

#define spanX (80.0*scaleFactor)
#define spanY (80.0*scaleFactor)
#define spanZ (80.0*scaleFactor)

// =================================================================================
// Set the element parameters
// =================================================================================
// The number of elements in each direction is 2^(refineFactor) * subdivisions
// For optimal performance, use refineFactor primarily to determine the element size
#define subdivisionsX 3 //10
#define subdivisionsY 3 //10
#define subdivisionsZ 1 //5
#define refineFactor 5 //6

// Set the polynomial degree of the element (suggested values: 1 or 2)
#define finiteElementDegree 2

// =================================================================================
// Set the adaptive mesh refinement parameters
// =================================================================================
// Set the flag determining if adaptive meshing is activated
#define hAdaptivity true

// Set the maximum and minimum level of refinement
#define maxRefinementLevel (refineFactor+3)
#define minRefinementLevel (2)

// Set the fields used to determine the refinement. Fields determined by the order
// declared in "equations.h", starting at zero
#define refineCriterionFields {1,2}

// Set the maximum and minimum value of the fields where the mesh should be refined
#define refineWindowMax {0.99,0.99}
#define refineWindowMin {0.001,0.001}

// Set the number of time steps between remeshing operations
#define skipRemeshingSteps 1000

// =================================================================================
// Set the time step parameters
// =================================================================================
// The size of the time step
#define timeStep (1.0e-4*scaleFactor*scaleFactor) //1.1e-4*scaleFactor*scaleFactor
#define timeIncrements 1000000 //5000000
#define timeFinal (timeStep*timeIncrements)


// =================================================================================
// Set the elliptic solver parameters
// =================================================================================
// The solver type (currently the only recommended option is conjugate gradient)
#define solverType "SolverCG"

// The flag that determines whether the tolerance for solver convergence should
// be an absolute tolerance (absTol=true) or a relative tolerance (absTol=false)
#define absTol true

// The tolerance for convergence (L2 norm of the residual)
#define solverTolerance 1.0e-1

// The maximum number of solver iterations per time step
#define maxSolverIterations 10000

// =================================================================================
// Set the output parameters
// =================================================================================
// Each field in the problem will be output is writeOutput is set to "true"
#define writeOutput true

// Type of spacing between outputs ("EQUAL_SPACING", "LOG_SPACING", "N_PER_DECADE",
// or "LIST")
#define outputCondition "EQUAL_SPACING"

// Number of times the program outputs the fields (total number for "EQUAL_SPACING"
// and "LOG_SPACING", number per decade for "N_PER_DECADE", ignored for "LIST")
#define numOutputs 20

// User-defined list of time steps where the program should output. Only used if
// outputCondition is "LIST"
#define outputList {0}

#define skipPrintSteps 1000

// =================================================================================
// Set the flag determining if the total free energy is calculated for each output
// =================================================================================
#define calcEnergy true

#define nucleationOccurs true









