//Parameters list for beta prime precipitation evolution problem

//define problem dimensions
#define problemDIM 2
#define spanX 40.0 //14.0
#define spanY 40.0 //14.0
#define spanZ 40.0 //10.0 //14.0

//define mesh parameters
#define subdivisionsX 3
#define subdivisionsY 3
#define subdivisionsZ 3
#define refineFactor 5
#define finiteElementDegree 2

// Define number of fields in the problem
// n1, n2, n3, c, u
// Cahn Hilliard part has no gradient term,
// hence chemical potential (mu) field not required as mixed formulation is not needed.
#define num_sop 3							// for now, must be between 1 and 3
#define numFields (1+num_sop+problemDIM)

//define time step parameters
#define timeStep (4.0e-4)
#define timeIncrements 5000
#define timeFinal 100.0
#define skipImplicitSolves 1

//define solver parameters
#define solverType SolverCG
#define abs_tol true
#define relSolverTolerance 1.0e-2
#define absSolverTolerance 1.0e-4
#define maxSolverIterations 1000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 100


