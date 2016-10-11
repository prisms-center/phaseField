//Parameter list for Mechanics (infinitesimal strain) problem

//define problem dimensions
#define problemDIM 3
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 4
#define finiteElementDegree 1

//define number of fields in the problem
//ux, uy, uz
#define numFields 3

//define solver paramters
#define solverType SolverCG
#define relSolverTolerance 1.0e-7
#define maxSolverIterations 10000

//define material properties 
#define MaterialModelV ISOTROPIC
#define MaterialConstantsV {1.0,0.3}

//adaptive refinement parameters
#define hAdaptivity true
#define maxRefinementLevel (refineFactor+1)
#define minRefinementLevel (0)
#define refinementDOF 0
#define topRefineFraction 0.99
#define bottomCoarsenFraction 0.01
