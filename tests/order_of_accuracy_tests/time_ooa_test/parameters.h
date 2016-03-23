//Parameters list for beta prime precipitation evolution problem

//define problem dimensions
#define problemDIM 3
#define spanX 32.0
#define spanY 32.0
#define spanZ 32.0

//define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 4
#define finiteElementDegree 1

//define time step parameters
//#define timeStep 3.0e-5
//#define timeIncrements 10000
//#define timeFinal 100000000
//#define skipImplicitSolves 10000000

//define solver parameters
#define solverType SolverCG
#define relSolverTolerance 1.0e-4
#define maxSolverIterations 10000

//define results output parameters
#define writeOutput true
#define skipOutputSteps (timeIncrements/10)

// flag to allow or disallow nucleation
#define nucleation_occurs false

#define numFields (4+problemDIM)

//define Cahn-Hilliard parameters (No Gradient energy)
#define McV 1.0

//define Allen-Cahn parameters
#define Mn1V 50.0
#define Mn2V 50.0
#define Mn3V 50.0

double Kn1[3][3]={{1.0,0,0},{0,0.5,0},{0,0,1.0}};
double Kn2[3][3]={{1.0,0,0},{0,0.5,0},{0,0,1.0}};
double Kn3[3][3]={{1.0,0,0},{0,0.5,0},{0,0,1.0}};

//define energy barrier coefficient (used to tune the interfacial energy)
#define W -1.0

//define Mechanical properties
#define MaterialModelV ISOTROPIC
#define MaterialConstantsV {1.0,0.3}
double sf1Strain[3][3] = {{0,0,0},{0,0,0},{0,0,0}}; //Mg-Nd beta-prime
double sf2Strain[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sf3Strain[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

// --------------------------------------------
// Unused for these tests
#define n_dependent_stiffness false
#define c_dependent_misfit false
// Used in n=1 phase if n_dependent_stiffness == true
#define MaterialModelBetaV ISOTROPIC
#define MaterialConstantsBetaV {2.0,0.3}

// Stress-free transformation strains (concentration dependent terms, used if c_dependent_misfit == true)
// Linear fits for the stress-free transformation strains in for sfts_p = ap * c + bp
double a1[3][3] = {{0.1,0,0},{0,0.6,0},{0,0,0}};
double b1[3][3] = {{-0.01,0,0},{0,-0.1,0},{0,0,0}};

double a2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double b2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

double a3[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double b3[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
// --------------------------------------------

//define free energy expressions (Mg-Nd data from CASM)
#define faV (24.7939*c*c - 1.6752*c - 1.9453e-06)
#define facV (49.5878*c - 1.6752)
#define faccV (49.5878)
#define fbV (37.9316*c*c - 10.7373*c + 0.5401)
#define fbcV (75.8633*c - 10.7373)
#define fbccV (75.8633)

#define h1V (3.0*n1*n1-2.0*n1*n1*n1)
#define h2V (3.0*n2*n2-2.0*n2*n2*n2)
#define h3V (3.0*n3*n3-2.0*n3*n3*n3)
#define hn1V (6.0*n1-6.0*n1*n1)
#define hn2V (6.0*n2-6.0*n2*n2)
#define hn3V (6.0*n3-6.0*n3*n3)

// This double-well function can be used to tune the interfacial energy
#define fbarrierV (n1*n1-2.0*n1*n1*n1+n1*n1*n1*n1)
#define fbarriernV (2.0*n1-6.0*n1*n1+4.0*n1*n1*n1)

// Residuals
#define rcV   (c)
#define rcxTemp ( cx*((1.0-h1V-h2V-h3V)*faccV+(h1V+h2V+h3V)*fbccV) + n1x*((fbcV-facV)*hn1V) + n2x*((fbcV-facV)*hn2V) + n3x*((fbcV-facV)*hn3V) )
#define rcxV  (constV(-timeStep*McV)*rcxTemp)

#define rn1V   (n1-constV(timeStep*Mn1V)*((fbV-faV)*hn1V+W*fbarriernV+nDependentMisfitAC1))
#define rn2V   (n2-constV(timeStep*Mn2V)*((fbV-faV)*hn2V))
#define rn3V   (n3-constV(timeStep*Mn3V)*((fbV-faV)*hn3V))
#define rn1xV  (constV(-timeStep*Mn1V)*Knx1)
#define rn2xV  (constV(-timeStep*Mn2V)*Knx2)
#define rn3xV  (constV(-timeStep*Mn3V)*Knx3)

// Initial geometry
#define x_denom 1.0
#define y_denom 1.0
#define z_denom 1.0
#define initial_interface_coeff 1.0
#define initial_radius 7.0
#define avg_Nd 0.004
