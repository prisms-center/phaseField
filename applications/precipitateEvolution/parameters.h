//Parameter list for the Beta Prime precipitate evolution problem 
//(Coupled Allen Cahn, Cahn Hilliard and Mechanics formulation)
//The free energy expressions in this file are from the reference:
//H. Liu et al, "A simulation study of the shape of beta prime precipitates in Mg–Y and Mg–Gd alloys", 
//Acta Materialia, Volume 61, Issue 2, January 2013, Pages 453-466. http://dx.doi.org/10.1016/j.actamat.2012.09.044

// Define problem dimensions
#define problemDIM 2
#define spanX 40.0
#define spanY 40.0
#define spanZ 40.0

// Define mesh parameters
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

// Define time step parameters
#define timeStep 4.0e-4
#define timeFinal 100.0
#define timeIncrements 5000
#define skipImplicitSolves 1

// Define solver paramters
#define solverType SolverCG
#define abs_tol true
#define relSolverTolerance 1.0e-2
#define absSolverTolerance 1.0e-4
#define maxSolverIterations 1000

// Define results output parameters
#define writeOutput true
#define skipOutputSteps 500

#define calc_energy true





