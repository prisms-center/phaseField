// Parameter list for the coupled Allen-Cahn/Cahn-Hilliard example application
// All strictly numerical parameters should be set in this file

// =================================================================================
// Set the number of dimensions (1, 2, or 3 for a 1D, 2D, or 3D calculation)
// =================================================================================
#define problemDIM 2

// =================================================================================
// Set the length of the domain in all three dimensions
// =================================================================================
// Each axes spans from zero to the specified length
#define spanX 800.0
#define spanY 800.0
#define spanZ 800.0

// =================================================================================
// Set the element parameters
// =================================================================================
// The number of elements in each direction is 2^(refineFactor) * subdivisions
// For optimal performance, use refineFactor primarily to determine the element size
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 10

// Set the polynomial degree of the element (suggested values: 1 or 2)
#define finiteElementDegree 2

// =================================================================================
// Set the time step parameters
// =================================================================================
// The size of the time step
#define timeStep 0.011

// The simulation ends when either timeFinal is reached or the number of time steps
// equals timeIncrements
#define timeFinal 440.0
#define timeIncrements 40000

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
#define calcEnergy true

// =================================================================================
// Set the KKS model parameters
// =================================================================================

//Average alloy composition (for initial contitions)
#define c_avg 0.2

// Cahn-Hilliard mobility
#define McV 1.0

// Allen-Cahn mobility
#define MnV 1.0

// Allen-Cahn gradient energy coefficient
#define KnV 2.0

// Energy barrier coefficient (used to tune the interfacial energy)
#define W (KnV)

// Interface coefficient
double interface_coeff=sqrt(2.0*KnV/W);

//Constants for free energy expressions
double A0 = 0.0;
double A2 = 1.0;
double calmin= 0.0;
double B0 = 0.0;
double B2 = 4.0;
double cbtmin = 1.0;

// =================================================================================
// Set the nucleation parameters
// =================================================================================

// Flag for nucleation to occur
#define nucleation_occurs true

// Nucleation radius (order parameter)
#define n_radius 5.0

// Hold time for order parameter
#define t_hold 20.0

// Small constant for sign function
#define epsil 1.0e-7

// Minimum distance between nuclei
#define minDistBetwenNuclei (4.0*n_radius)

// Number of time steps between nucleation attempts
#define skipNucleationSteps 30

// radius for order parameter hold
#define opfreeze_radius (2.0*n_radius)

//Minimum distance from the edges of the system where nucleation can occur
#define borderreg (2.0*n_radius)

// Constants k1 and k2 for nucleation rate in the bulk
#define k1 498.866
#define k2 4.14465




