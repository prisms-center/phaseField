//Parameter list for Mechanics (infinitesimal strain) problem

//define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh parameters
#define refineFactor 6
#define finiteElementDegree 1

//define time step parameters
#define dt 1.0e-2
#define numIncrements 10

//define results output parameters
#define writeOutput true
#define skipOutputSteps 1

//define material properties
#define Ev 1.0
#define nuV 0.3
#define muV (Ev/(2*(1+nuV)))
#define lambdaV (nuV*Ev/((1+nuV)*(1-2*nuV))) 
#define CijklV (muV*((i==k)*(j==l)+(i==l)*(j==k))+lambdaV*(i==j)*(k==l))
//define required residuals
#define ruxV  (ux)

