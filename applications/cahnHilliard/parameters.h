//Parameter list for the Cahn-Hilliard spinodal decomposition problem

//define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh parameters
#define refineFactor 7
#define finiteElementDegree 1

//define time step parameters
#define dt 1.0e-2
#define numIncrements 1000000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 100000

//define Cahn-Hilliard parameters
#define Mc 1.0 
#define Kc 1.0 

//define free energy expressions
#define fcV (4.0*c*(c-1.0)*(c-0.5))

//define required residuals
#define rmuV  (fcV)
#define rmuxV (cx)
#define rcV   (c)
#define rcxV  (mux)

