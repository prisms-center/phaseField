//Parameter list for the Allen-Cahn problem

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
//n
#define numFields 1

//define time step parameters
#define timeStep 1.0e-3
#define timeFinal 20.0
#define timeIncrements 200000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 1000


