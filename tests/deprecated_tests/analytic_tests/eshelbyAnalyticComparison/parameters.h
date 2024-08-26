// Parameter list for the Beta Prime precipitate evolution problem
//(Coupled Allen Cahn, Cahn Hilliard and Mechanics formulation)
// The free energy expressions in this file are from the reference:
// H. Liu et al, "A simulation study of the shape of beta prime precipitates in
// Mg–Y and Mg–Gd alloys", Acta Materialia, Volume 61, Issue 2, January 2013,
// Pages 453-466. http://dx.doi.org/10.1016/j.actamat.2012.09.044

// Define problem dimensions
#define problemDIM 3
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

// Define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 6
#define finiteElementDegree 1

// define number of fields in the problem
#define numFields 3

// Define solver paramters
#define solverType SolverCG
#define absTol false
#define solverTolerance 1.0e-10
#define maxSolverIterations 1000

// Define results output parameters
#define writeOutput true

#define calcEnergy false
