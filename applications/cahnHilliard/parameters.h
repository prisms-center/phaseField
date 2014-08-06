//Parameter list for the Cahn-Hilliard spinodal decomposition problem

//define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh paramters
#define refineFactor 7
#define finiteElementDegree 1

//define time step paramters
#define dt 5.0e-3
#define numIncrements 1000

//define results output paramters
#define writeOutput true
#define skipOutputSteps 100

//define Cahn-Hilliard paramters
#define Mc 1.0 
#define Kc 1.0 

//define free energy expressions
#define fcV (4.0*c*(c-1.0)*(c-0.5))

//define required residuals
#define rmuV  (fcV)
#define rmuxV (cx)
#define rcV   (c)
#define rcxV  (mux)

