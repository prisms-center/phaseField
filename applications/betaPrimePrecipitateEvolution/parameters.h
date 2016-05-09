//Parameter list for the Beta Prime precipitate evolution problem 
//(Coupled Allen Cahn, Cahn Hilliard and Mechanics formulation)
//The free energy expressions in this file are from the reference:
//H. Liu et al, "A simulation study of the shape of beta prime precipitates in Mg–Y and Mg–Gd alloys", 
//Acta Materialia, Volume 61, Issue 2, January 2013, Pages 453-466. http://dx.doi.org/10.1016/j.actamat.2012.09.044

// Define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

// Define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 7
#define finiteElementDegree 1

// Define number of fields in the problem
// n1, n2, n3, c, u
// Cahn Hilliard part has no gradient term,
// hence chemical potential (mu) field not required as mixed formulation is not needed.
#define num_sop 3							// for now, must be between 1 and 3
#define numFields (1+num_sop+problemDIM)

// flag to allow or disallow nucleation
#define nucleation_occurs false

// Define time step parameters
#define timeStep 1.0e-4
#define timeFinal 10.0
#define timeIncrements 10000
#define skipImplicitSolves 1 //1000

// Define solver paramters
#define solverType SolverCG
#define relSolverTolerance 1.0e-4
#define maxSolverIterations 1000

// Define results output parameters
#define writeOutput true
#define skipOutputSteps 1000





