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
#define relSolverTolerance 1.0e-10
#define maxSolverIterations 1000

//define material properties 
#define MaterialModelV ISOTROPIC
#define MaterialConstantsV {1.0,0.3}
