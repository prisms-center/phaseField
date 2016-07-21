//Parameter list for the Beta Prime precipitate evolution problem 
//(Coupled Allen Cahn, Cahn Hilliard and Mechanics formulation)
//The free energy expressions in this file are from the reference:
//H. Liu et al, "A simulation study of the shape of beta prime precipitates in Mg–Y and Mg–Gd alloys", 
//Acta Materialia, Volume 61, Issue 2, January 2013, Pages 453-466. http://dx.doi.org/10.1016/j.actamat.2012.09.044

// Define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

// Define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 7
#define finiteElementDegree 1

//define number of fields in the problem
//n
#define numFields 1

//define time step parameters
#define timeStep 1.0e-3
#define timeFinal 20.0
#define timeIncrements 20000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 1000

#define calc_energy true


// Definition of the variables in the model
#define num_var 1
std::string var_name[num_var] = {"n"};
std::string var_type[num_var] = {"SCALAR"};
std::string var_eq_type[num_var] = {"PARABOLIC"};
bool need_value[num_var] = {true};
bool need_gradient[num_var] = {true};
bool need_hessian[num_var] = {false};
bool value_residual[num_var] = {true};
bool gradient_residual[num_var] = {true};






