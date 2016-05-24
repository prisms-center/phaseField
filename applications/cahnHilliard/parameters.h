//Parameter list for the Cahn-Hilliard spinodal decomposition problem

//define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1

//define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 6
#define finiteElementDegree 1

//adaptive refinement parameters
#define hAdaptivity true
#define maxRefinementLevel refineFactor+1
#define refinementDOF 1 //refining using composition field
#define topRefineFraction 0.8
#define bottomCoarsenFraction 0.1


//define number of fields in the problem
//c, mu
#define numFields 2

//define time step parameters
#define timeStep 1.0e-3
#define timeFinal 1000
#define timeIncrements 1000000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 10000

//define Cahn-Hilliard parameters
#define McV 1.0 
#define KcV 4.0 

//define free energy expressions
#define fcV (4.0*c*(c-1.0)*(c-0.5))

//define required residuals
#define rmuV  (fcV)
#define rmuxV (constV(KcV)*cx)
#define rcV   (c)
#define rcxV  (constV(-McV*timeStep)*mux)


