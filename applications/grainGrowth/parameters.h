// Parameter list for the grain growth example problem.
// The model used for this example application is a simplified version of the model from:
// Gentry and Thornton, "Simulating Recrystalization in Titanium using the Phase Field Method",
// IOP Conference Series: Mater. Sci. Eng. 89 (2015) 012024.

// Define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

// Define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 7
#define finiteElementDegree 1

//define number of fields in the problem
//n
#define numFields 5

//define time step parameters
#define timeStep 1.0e-3
#define timeFinal 20.0
#define timeIncrements 20000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 1000

#define calcEnergy false







