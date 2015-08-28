//Parameter list for the Fickian diffusion problem

//define problem dimensions
#define problemDIM 2
#define spanX 1.0
#define spanY 1.0
#define spanZ 1.0

//define mesh parameters
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

//define Fickian diffusion parameters
#define DcV 0.01

//define geometric expression
//defining a circle centered at (0.0,0.0) with radius spanX/10
#define exp1 (std::sqrt(x*x+y*y)<spanX/10)

//define temporal expression
#define exp2 (t<timeFinal/4)


//define required residuals
#define rcV   (c + timeStep*subdomain(exp1 && exp2, constV(1.0)) )
#define rcxV  (constV(-DcV*timeStep)*cx)


