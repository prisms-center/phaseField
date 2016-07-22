//Parameter list for the Fickian diffusion problem

// Define problem dimensions
#define problemDIM 2
#define spanX 1.0
#define spanY 1.0
#define spanZ 1.0

// Define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 7
#define finiteElementDegree 1

//define number of fields in the problem
//c
#define numFields 1

//define time step parameters
#define timeStep 1.0e-3
#define timeFinal 10.0
#define timeIncrements 10000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 100

#define calc_energy false


// Definition of the variables in the model
#define num_var 1
#define variable_name {"c"}
#define variable_type {"SCALAR"}
#define variable_eq_type {"PARABOLIC"}
#define need_val {true}
#define need_grad {true}
#define need_hess  {false}
#define need_val_residual {true}
#define need_grad_residual {true}






