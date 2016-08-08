// Parameter list for the precipitate evolution example application
// All strictly numerical parameters should be set in this file

// =================================================================================
// Set the number of dimensions (1, 2, or 3 for a 1D, 2D, or 3D calculation)
// =================================================================================
#define problemDIM 2
// =================================================================================

// =================================================================================
// Set the length of the domain in all three dimensions
// =================================================================================
// Each axes spans from zero to the specified length
#define spanX 40.0
#define spanY 40.0
#define spanZ 40.0
// =================================================================================

// =================================================================================
// Set the element parameters
// =================================================================================
// The number of elements in each direction is 2^(refineFactor) * subdivisions
// For optimal performance, use refineFactor primarily to determine the element size
#define subdivisionsX 3
#define subdivisionsY 3
#define subdivisionsZ 3
#define refineFactor 5

// Set the polynomial degree of the element
// Suggested values are either 1 or 2
#define finiteElementDegree 2
// =================================================================================

// =================================================================================
// Set the time step parameters
// =================================================================================
// The size of the time step
#define timeStep 4.0e-4
#define timeFinal 100.0
#define timeIncrements 5000
// =================================================================================

// =================================================================================
// Set the elliptic solver parameters
// =================================================================================
// Set the solver type (currently the only recommended option is conjugate gradient)
#define solverType SolverCG

// Set the flag that determines whether the tolerance for solver convergence should
// be an absolute tolerance (abs_tol=true) or a relative tolerance (abs_tol=false)
#define abs_tol true

// Relative and absolute tolerances for convergence (L2 norm of the residual)
#define relSolverTolerance 1.0e-2
#define absSolverTolerance 1.0e-4

// Maximum number of solver iterations per time step
#define maxSolverIterations 1000
// =================================================================================

// =================================================================================
// Set the output parameters
// =================================================================================
// Each field in the problem will be output is writeOutput is set to "true"
#define writeOutput true

// Output files are written every skipOutputSteps time steps
#define skipOutputSteps 500
// =================================================================================

// =================================================================================
// Set the flag determining if the total free energy is calculated for each output
// =================================================================================
#define calc_energy true
// =================================================================================





