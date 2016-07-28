//Parameter list for the coupled Cahn-Hilliard and Mechanics problem
// General interface
// Define problem dimensions
#define problemDIM 2
#define spanX 400.0
#define spanY 100.0
#define spanZ 40.0

// Define mesh parameters
#define subdivisionsX 4
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 6
#define finiteElementDegree 2

// Define time step parameters
#define timeStep 1.0e-3
#define timeFinal 5000.0
#define timeIncrements 5e6
#define skipImplicitSolves 1

// Define solver parameters
#define solverType SolverCG
#define relSolverTolerance 1.0e-5
#define maxSolverIterations 1000
#define abs_tol true
#define absSolverTolerance 1.0e-3

// Define results output parameters
#define writeOutput true
#define skipOutputSteps 1000

// Definition of the variables in the model
#define num_var 3 // total number of variable 
#define variable_name {"c", "mu", "u"}
#define variable_type {"SCALAR","SCALAR","VECTOR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC","ELLIPTIC"}
#define need_val {true, false, false}
#define need_grad {true, true, true}
#define need_hess {false, false, false} // Currently overridden based on value of "n_dependent_stiffness"
#define need_val_residual {true, true, false}
#define need_grad_residual {true, true, true} // ith term for ith Equ

// LHS denotes the CG-solver 
#define need_val_LHS {true, false, false}
#define need_grad_LHS {false, false, true}
#define need_hess_LHS {false, false, false}
#define need_val_residual_LHS {false, false, false}
#define need_grad_residual_LHS {false, false, true}

// Define Cahn-Hilliard parameters (no gradient energy terms)
#define McV 1.0 
#define KcV 0.1  

// Define Mechanical properties
#define c_dependent_stiffness true
#define n_dependent_stiffness true  // DELETE 
// Mechanical symmetry of the material and stiffness parameters
// Used throughout system if c_dependent_stiffness == false, used in c=0 phase if c_dependent_stiffness == true
#define MaterialModel {{ANISOTROPIC},{ANISOTROPIC}}
// #define MaterialConstantsV {1.0,1.0,0.385,0.5,0.0,0.0}
 #define MaterialConstants {{10.0,1.0,0.385,1.0,0.0,0.0},{10.0,1.0,0.385,1.0,0.0,0.0}}

//2D models:
//ISOTROPIC- (Plane Strain) 2 constants [E, nu]
//ANISOTROPIC- 6 constants [C11 C22 C33 C12 C13 C23]

// Used in n=1 phase if n_dependent_stiffness == true
//#define MaterialModelBetaV ANISOTROPIC
// #define MaterialConstantsBetaV {1.0,1.0,0.385,0.5,0.0,0.0}
//#define MaterialConstantsBetaV {10.0,1.0,0.385,1.0,0.0,0.0}

// Stress-free transformation strains
// Linear fits for the stress-free transformation strains in for sfts = sfts_linear * c + sfts_const
double sfts_const1[3][3] = {{0.00,0.00,0},{0,0.0,0},{0,0,0.00}};
double sfts_const2[3][3] = {{0.05,0.00,0},{0,0.01,0},{0,0,0.00}};

//define free energy expressions  
#define fV (0.25*c*c*(1.0-c)*(1.0-c)) // THIS IS USED JUST TO REPRESENT THE residuals below in a short form
#define fcV (0.5*(c*(1.0-c)*(1.0-c) -c*c*(1.0-c)))  // THIS IS USED JUST TO REPRESENT THE residuals below  in a short form

//define required residuals  
#define rmuV  (fcV -cDependentMisfitCH)
#define rmuxV (constV(KcV)*cx)
#define rcV   (c)
#define rcxV  (constV(-McV*timeStep)*mux)

#define calc_energy  true
