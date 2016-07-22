//Parameter list for the Beta Prime precipitate evolution problem 
//(Coupled Allen Cahn, Cahn Hilliard and Mechanics formulation)
//The free energy expressions in this file are from the reference:
//H. Liu et al, "A simulation study of the shape of beta prime precipitates in Mg–Y and Mg–Gd alloys", 
//Acta Materialia, Volume 61, Issue 2, January 2013, Pages 453-466. http://dx.doi.org/10.1016/j.actamat.2012.09.044

// Define problem dimensions
#define problemDIM 2
#define spanX 40.0
#define spanY 40.0
#define spanZ 40.0

// Define mesh parameters
#define subdivisionsX 3
#define subdivisionsY 3
#define subdivisionsZ 3
#define refineFactor 5
#define finiteElementDegree 2

// Define number of fields in the problem
// n1, n2, n3, c, u
// Cahn Hilliard part has no gradient term,
// hence chemical potential (mu) field not required as mixed formulation is not needed.
#define num_sop 3							// for now, must be between 1 and 3
#define numFields (1+num_sop+problemDIM)

// Define time step parameters
#define timeStep 4.0e-4
#define timeFinal 100.0
#define timeIncrements 5000
#define skipImplicitSolves 1

// Define solver paramters
#define solverType SolverCG
#define abs_tol true
#define relSolverTolerance 1.0e-2
#define absSolverTolerance 1.0e-4
#define maxSolverIterations 1000

// Define results output parameters
#define writeOutput true
#define skipOutputSteps 500

// Definition of the variables in the model
#define num_var 5
std::string var_name[num_var] = {"c", "n1", "n2", "n3", "u"};
std::string var_type[num_var] = {"SCALAR","SCALAR","SCALAR","SCALAR","VECTOR"};
std::string var_eq_type[num_var] = {"PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","ELLIPTIC"};
bool need_value[num_var] = {true, true, true, true, false};
bool need_gradient[num_var] = {true, true, true, true, true};
bool need_hessian[num_var] = {false, false, false, false, false}; // Currently overridden based on value of "n_dependent_stiffness"
bool value_residual[num_var] = {true, true, true, true, false};
bool gradient_residual[num_var] = {true, true, true, true, true};

//bool need_value_LHS[num_var] = {false, true, true, true, false};
//bool need_gradient_LHS[num_var] = {false, false, false, false, true};
//bool need_hessian_LHS[num_var] = {false, false, false, false, false};
//bool value_residual_LHS[num_var] = {false, false, false, false, false};
//bool gradient_residual_LHS[num_var] = {false, false, false, false, true};

#define need_val_LHS {false, true, true, true, false}
#define need_grad_LHS {false, false, false, false, true}
#define need_hess_LHS {false, false, false, false, false}
#define need_val_residual_LHS {false, false, false, false, false}
#define need_grad_residual_LHS {false, false, false, false, true}

// Define Cahn-Hilliard parameters (no gradient energy terms)
#define McV 1.0

// Define Allen-Cahn parameters
#define Mn1V 100.0
#define Mn2V 100.0
#define Mn3V 100.0

double Kn1[3][3]={{0.03,0,0},{0,0.007,0},{0,0,1.0}};
double Kn2[3][3]={{0.01275,-0.009959,0},{-0.009959,0.02425,0},{0,0,1.0}};
double Kn3[3][3]={{0.01275,0.009959,0},{0.009959,0.02425,0},{0,0,1.0}};

// Define Mechanical properties
#define n_dependent_stiffness false
// Mechanical symmetry of the material and stiffness parameters
// Used throughout system if n_dependent_stiffness == false, used in n=0 phase if n_dependent_stiffness == true
#define MaterialModelV ISOTROPIC
#define MaterialConstantsV {2.0,0.3}

// Used in n=1 phase if n_dependent_stiffness == true
#define MaterialModelBetaV ISOTROPIC
#define MaterialConstantsBetaV {2.0,0.3}

// Stress-free transformation strains
// Linear fits for the stress-free transformation strains in for sfts = sfts_linear * c + sfts_const
double sfts_linear1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sfts_const1[3][3] = {{0.0345,0,0},{0,0.0185,0},{0,0,-0.00270}};

double sfts_linear2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sfts_const2[3][3]={{0.0225,-0.0069,0},{-0.0069,0.0305,0},{0,0,-0.00270}};

double sfts_linear3[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sfts_const3[3][3]={{0.0225, 0.0069,0},{0.0069,0.0305,0},{0,0,-0.00270}};

//define free energy expressions
#define faV (-1.6704-4.776*c+5.1622*c*c-2.7375*c*c*c+1.3687*c*c*c*c)
#define facV (-4.776 + 10.3244*c - 8.2125*c*c + 5.4748*c*c*c)
#define faccV (10.3244-16.425*c+16.4244*c*c)
#define fbV (5.0*c*c-5.9746*c-1.5924)
#define fbcV (10.0*c-5.9746)
#define fbccV (10.0)
#define h1V (10.0*n1*n1*n1-15.0*n1*n1*n1*n1+6.0*n1*n1*n1*n1*n1)
#define h2V (10.0*n2*n2*n2-15.0*n2*n2*n2*n2+6.0*n2*n2*n2*n2*n2)
#define h3V (10.0*n3*n3*n3-15.0*n3*n3*n3*n3+6.0*n3*n3*n3*n3*n3)
#define hn1V (30.0*n1*n1-60.0*n1*n1*n1+30.0*n1*n1*n1*n1)
#define hn2V (30.0*n2*n2-60.0*n2*n2*n2+30.0*n2*n2*n2*n2)
#define hn3V (30.0*n3*n3-60.0*n3*n3*n3+30.0*n3*n3*n3*n3)

#define calc_energy true





