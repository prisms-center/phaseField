//Parameter list for the Allen-Cahn problem

//define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh parameters
#define refineFactor 7
#define finiteElementDegree 1

//define number of fields in the problem
//n
#define numFieldsV 1

//define time step parameters
#define timeStepV 1.0e-3
#define finalTimeV 20.0
#define totalIncrementsV 200000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 1000

//define Allen-Cahn parameters
#define MnV 1.0 
#define KnV 4.0 

//define free energy expressions
#define fnV (4.0*n*(n-1.0)*(n-0.5))

//define required residuals
#define rnV  (n-(dt*MnV)*fnV)
#define rnxV (nx)

