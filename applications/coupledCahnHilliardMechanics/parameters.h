//Parameter list for the coupled Cahn-Hilliard and Mechanics problem

//define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh parameters
#define refineFactor 7
#define finiteElementDegree 1

//define number of fields in the problem
//c, mu, u
#define numFields (2+problemDIM)

//define time step parameters
#define timeStep 1.0e-2
#define timeFinal 100.0
#define timeIncrements 10000
#define skipImplicitSolves 1000

//define solver paramters
#define solverType SolverCG
#define relSolverTolerance 1.0e-10
#define maxSolverIterations 1000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 1000

//define Cahn-Hilliard parameters
#define McV 1.0 
#define KcV 1.0 

//define material properties 
#define MaterialModelV ISOTROPIC
#define MaterialConstantsV {1.0,0.3}

//define free energy expressions
#define fcV (4.0*c*(c-1.0)*(c-0.5))

//define required residuals
#define rmuV  (fcV)
#define rmuxV (constV(KcV)*cx)
#define rcV   (c)
#define rcxV  (constV(-McV*timeStep)*mux)

//define chemical strain
#define chemicalStrainV (-0.01*c)


