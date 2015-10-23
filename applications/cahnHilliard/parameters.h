//Parameter list for the Cahn-Hilliard spinodal decomposition problem

//define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 7
#define finiteElementDegree 1

//define number of fields in the problem
//c, mu
#define numFields 2

//define time step parameters
#define timeStep 1.0e-2
#define timeFinal 100.0
#define timeIncrements 10000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 100

//define Cahn-Hilliard parameters
#define McV 1.0 
#define KcV 1.0 

//define free energy expressions
#define fcV (4.0*c*(c-1.0)*(c-0.5))

//define required residuals
#define rmuV  (fcV)
#define rmuxV (constV(KcV)*cx)
#define rcV   (c)
#define rcxV  (constV(-McV*timeStep)*mux)


