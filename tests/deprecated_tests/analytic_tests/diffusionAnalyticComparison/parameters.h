// Parameter list for the Fickian diffusion problem

// Define problem dimensions
#define problemDIM 2
#define spanX 10.0
#define spanY 1.0
#define spanZ 10.0

// Define mesh parameters
#define subdivisionsX 500
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 0
#define finiteElementDegree 1

// define number of fields in the problem
// c
#define numFields 1

// define time step parameters
#define timeStep 1.0e-4
#define timeFinal 100.0
#define timeIncrements 10000

// define results output parameters
#define writeOutput true

// Type of spacing between outputs ("EQUAL_SPACING", "LOG_SPACING", or
// "N_PER_DECADE")
#define outputCondition "EQUAL_SPACING"

// Number of times the program outputs the fields (total number for
// "EQUAL_SPACING" and "LOG_SPACING", number per decade for "N_PER_DECADE")
#define numOutputs 10

#define calcEnergy false
