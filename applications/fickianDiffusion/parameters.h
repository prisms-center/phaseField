// Parameter list for the diffusion example application
// All strictly numerical parameters should be set in this file

// =================================================================================
// Set the number of dimensions (1, 2, or 3 for a 1D, 2D, or 3D calculation)
// =================================================================================
#define problemDIM 2

// =================================================================================
// Set the length of the domain in all three dimensions
// =================================================================================
// Each axes spans from zero to the specified length
#define spanX 1.0
#define spanY 1.0
#define spanZ 1.0

// =================================================================================
// Set the element parameters
// =================================================================================
// The number of elements in each direction is 2^(refineFactor) * subdivisions
// For optimal performance, use refineFactor primarily to determine the element size
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 7

// Set the polynomial degree of the element (suggested values: 1 or 2)
#define finiteElementDegree 1

// =================================================================================
// Set the time step parameters
// =================================================================================
// The size of the time step
#define timeStep 1.0e-3

// The simulation ends when either timeFinal is reached or the number of time steps
// equals timeIncrements
#define timeFinal 10.0
#define timeIncrements 10000

// =================================================================================
// Set the output parameters
// =================================================================================
// Each field in the problem will be output is writeOutput is set to "true"
#define writeOutput true

// Output files are written every skipOutputSteps time steps
#define skipOutputSteps 100

// =================================================================================
// Set the flag determining if the total free energy is calculated for each output
// =================================================================================
#define calcEnergy false








