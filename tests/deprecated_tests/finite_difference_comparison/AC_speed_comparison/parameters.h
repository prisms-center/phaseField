// Parameter list for the Beta Prime precipitate evolution problem
//(Coupled Allen Cahn, Cahn Hilliard and Mechanics formulation)
// The free energy expressions in this file are from the reference:
// H. Liu et al, "A simulation study of the shape of beta prime precipitates in
// Mg–Y and Mg–Gd alloys", Acta Materialia, Volume 61, Issue 2, January 2013,
// Pages 453-466. http://dx.doi.org/10.1016/j.actamat.2012.09.044

// Define problem dimensions
#define problemDIM 3
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

// Define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 7
#define finiteElementDegree 1

// define number of fields in the problem
// n
#define numFields 1

// define time step parameters
#define timeStep 1.0e-3
#define timeFinal 20.0
#define timeIncrements 30000

// define results output parameters
#define writeOutput true

// Type of spacing between outputs ("EQUAL_SPACING", "LOG_SPACING", or
// "N_PER_DECADE")
#define outputCondition "EQUAL_SPACING"

// Number of times the program outputs the fields (total number for
// "EQUAL_SPACING" and "LOG_SPACING", number per decade for "N_PER_DECADE")
#define numOutputs 1

#define calcEnergy false
