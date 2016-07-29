//Parameter list for the coupled Cahn-Hilliard and Mechanics problem
// General interface
// Define problem dimensions
#define problemDIM 2
#define spanX 400.0
#define spanY 100.0
#define spanZ 40.0

// Define mesh parameters
#define subdivisionsX 4
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 6
#define finiteElementDegree 2

// Define time step parameters
#define timeStep 1.0e-3
#define timeFinal 5000.0
#define timeIncrements 5e6
#define skipImplicitSolves 1

// Define solver parameters
#define solverType SolverCG
#define relSolverTolerance 1.0e-5
#define maxSolverIterations 1000
#define abs_tol true
#define absSolverTolerance 1.0e-3

// Define results output parameters
#define writeOutput true
#define skipOutputSteps 1000


#define calc_energy  false
