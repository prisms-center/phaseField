// Parameter list for the coupled Allen-Cahn/Cahn-Hilliard example application
// All strictly numerical parameters should be set in this file

// =================================================================================
// Set the number of dimensions (1, 2, or 3 for a 1D, 2D, or 3D calculation)
// =================================================================================
#define problemDIM 3
// =================================================================================

// =================================================================================
// Set the length of the domain in all three dimensions
// =================================================================================
// Each axes spans from zero to the specified length
#define spanX 64.0
#define spanY 64.0
#define spanZ 64.0
// =================================================================================

// =================================================================================
// Set the element parameters
// =================================================================================
// The number of elements in each direction is 2^(refineFactor) * subdivisions
// For optimal performance, use refineFactor primarily to determine the element
// size
#define subdivisionsX 5
#define subdivisionsY 5
#define subdivisionsZ 5
#define refineFactor 3

// Set the polynomial degree of the element
// Suggested values are either 1 or 2
#define finiteElementDegree 2
// =================================================================================

// =================================================================================
// Set the time step parameters
// =================================================================================
// The size of the time step
#define timeStep (5.0 / 5500.0)

// The simulation ends when either timeFinal is reached or the number of time
// steps equals timeIncrements
#define timeFinal 10000.0
#define timeIncrements 27500
// =================================================================================

// =================================================================================
// Set the output parameters
// =================================================================================
// Each field in the problem will be output is writeOutput is set to "true"
#define writeOutput true

// Type of spacing between outputs ("EQUAL_SPACING", "LOG_SPACING", or
// "N_PER_DECADE")
#define outputCondition "EQUAL_SPACING"

// Number of times the program outputs the fields (total number for
// "EQUAL_SPACING" and "LOG_SPACING", number per decade for "N_PER_DECADE")
#define numOutputs 1
// =================================================================================

// =================================================================================
// Set the flag determining if the total free energy is calculated for each
// output
// =================================================================================
#define calcEnergy false
// =================================================================================
