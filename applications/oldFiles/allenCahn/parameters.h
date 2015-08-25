//Parameter list for the Allen-Cahn problem

//define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh parameters
#define refineFactor 7
#define finiteElementDegree 1

//define time step parameters
#define dt 1.0e-4
#define numIncrements 200000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 1000

//define Allen-Cahn parameters
#define Mn 1.0 
#define Kn 4.0 

//define free energy expressions
#define fnV (4.0*n*(n-1.0)*(n-0.5))

//define required residuals
#define rnV  (n-(dt*Mn)*fnV)
#define rnxV (nx)

