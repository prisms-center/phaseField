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
#define refineFactor 5
#define finiteElementDegree 2

// Define time step parameters
#define timeStep 5e-2
#define timeFinal 5000.0
#define timeIncrements 100000 //5e6
#define skipImplicitSolves 1

// Define solver parameters
#define solverType SolverCG
#define maxSolverIterations 1000
#define absTol true
#define solverTolerance 1.0e-3

// Define results output parameters
#define writeOutput true

// Type of spacing between outputs ("EQUAL_SPACING", "LOG_SPACING", or "N_PER_DECADE")
#define outputCondition "EQUAL_SPACING"

// Number of times the program outputs the fields (total number for "EQUAL_SPACING"
// and "LOG_SPACING", number per decade for "N_PER_DECADE")
#define numOutputs 10


#define calc_energy  true
