//Parameter list for Mechanics (infinitesimal strain) problem

//define problem dimensions
#define problemDIM 3
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh parameters
#define refineFactor 4
#define finiteElementDegree 1

//define time step parameters
#define numIncrements 1

//define results output parameters
#define writeOutput true
#define skipOutputSteps 1

//define solver paramters
#define solverType SolverCG
#define relSolverTolerance 1.0e-10
#define maxSolverIterations 1000

//define material properties 
#define MaterialModelv ISOTROPIC
#define MaterialConstantsv {1.0,0.3}


