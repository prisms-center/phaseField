//Parameter list for the Beta Prime precipitate evolution problem 
//(Coupled Allen Cahn, Cahn Hilliard and Mechanics formulation)
//The free energy expressions in this file are from the reference:
//H. Liu et al, "A simulation study of the shape of beta prime precipitates in Mg–Y and Mg–Gd alloys", 
//Acta Materialia, Volume 61, Issue 2, January 2013, Pages 453-466. http://dx.doi.org/10.1016/j.actamat.2012.09.044

// Define problem dimensions
#define problemDIM 3
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

// Define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 4
#define finiteElementDegree 1

//define number of fields in the problem
#define numFields 3

// Define solver paramters
#define solverType SolverCG
#define relSolverTolerance 1.0e-10
#define maxSolverIterations 1000

// Define results output parameters
#define writeOutput true

// Definition of the variables in the model
#define num_var 1
#define variable_name {"u"}
#define variable_type {"VECTOR"}
#define variable_eq_type {"ELLIPTIC"}
#define need_val {false}
#define need_grad {true}
#define need_hess {false}
#define need_val_residual {false}
#define need_grad_residual {true}

#define need_val_LHS {false}
#define need_grad_LHS {true}
#define need_hess_LHS {false}
#define need_val_residual_LHS {false}
#define need_grad_residual_LHS {true}


// Define Mechanical properties
// Mechanical symmetry of the material and stiffness parameters
#define MaterialModelV ISOTROPIC
#define MaterialConstantsV {2.0,0.3}






