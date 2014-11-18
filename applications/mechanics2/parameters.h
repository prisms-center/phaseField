//Parameter list for the Cahn-Hilliard spinodal decomposition problem

//define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh parameters
#define refineFactor 2
#define finiteElementDegree 1

//define time step parameters
#define timeStepV 1.0e-2
#define finalTimeV 1.0
#define totalIncrementsV 10
 
//define results output parameters
#define writeOutput true
#define skipOutputSteps 1000

//define solver paramters
#define solverType SolverCG
#define relSolverTolerance 1.0e-10
#define maxSolverIterations 1000

//define material properties 
#define MaterialModelv ISOTROPIC
#define MaterialConstantsv {1.0,0.3}
