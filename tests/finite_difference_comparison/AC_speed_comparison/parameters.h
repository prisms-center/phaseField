//Parameter list for the Allen-Cahn problem

//define problem dimensions
#define problemDIM 3
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

// flag to allow or disallow nucleation
#define nucleation_occurs false

//define time step parameters
#define timeStep 1.0e-3
#define timeFinal 20.0
#define timeIncrements 30000

//define results output parameters
#define writeOutput true
#define skipOutputSteps (timeIncrements)

//define Allen-Cahn parameters
#define MnV 1.0 
#define KnV 4.0 

//define free energy expressions
#define fnV (4.0*n*(n-1.0)*(n-0.5))

//define required residuals
#define rnV  (n-constV(timeStep*MnV)*fnV)
#define rnxV (constV(-timeStep*KnV*MnV)*nx)

//define initial conditions
#define initial_radius 40.0
#define initial_interface_coeff 2.0

