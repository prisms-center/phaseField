//Parameter list for Mechanics (infinitesimal strain) problem

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
#define numIncrements 100

//define results output parameters
#define writeOutput true
#define skipOutputSteps 1

//define required residuals
#define ruxV  (ux)

