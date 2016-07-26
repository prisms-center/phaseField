//Parameter list for the Fickian diffusion problem

// Define problem dimensions
#define problemDIM 2
#define spanX 1.0
#define spanY 1.0
#define spanZ 1.0

// Define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 7
#define finiteElementDegree 1

//define number of fields in the problem
//c
#define numFields 1

//define time step parameters
#define timeStep 1.0e-3
#define timeFinal 10.0
#define timeIncrements 10000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 100

#define calc_energy false








