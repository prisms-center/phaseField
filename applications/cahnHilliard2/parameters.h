//Parameter list for the Cahn-Hilliard spinodal decomposition problem

//define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh parameters
#define refineFactor 7
#define finiteElementDegree 1

//define number of fields in the problem
//c, mu
#define numFieldsV 2

//define time step parameters
#define timeStepV 1.0e-2
#define finalTimeV 1000.0
#define totalIncrementsV 1000000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 1000

//define Cahn-Hilliard parameters
#define McV 1.0 
#define KcV 1.0 

//define free energy expressions
#define fcV (4.0*c*(c-1.0)*(c-0.5))

//define required residuals
#define rmuV  (fcV)
#define rmuxV (cx)
#define rcV   (c)
#define rcxV  (mux)

