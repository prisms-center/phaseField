//Parameters list for beta prime precipitation evolution problem

//define problem dimensions
#define problemDIM 2
#define spanX 43.0 //14.0
#define spanY 43.0 //14.0
#define spanZ 43.0 //10.0 //14.0

//define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 6
#define finiteElementDegree 2

//define time step parameters
#define timeStep (1.35e-5) //5e-6 //1.67e-5
#define timeIncrements 300000 //200000
#define timeFinal 100000000 //(timeStep*timeIncrements)
#define skipImplicitSolves 1

//define solver parameters
#define solverType SolverCG
#define abs_tol true
#define relSolverTolerance 1.0e-2
#define absSolverTolerance 1.0e-3
#define maxSolverIterations 10000

//define results output parameters
#define writeOutput true
#define skipOutputSteps (timeIncrements/10) //50000 //timeIncrements/10 //5000


