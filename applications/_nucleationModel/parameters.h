//Parameters list for the coupled Cahn Hilliard and Allen Cahn problem

//define problem dimensions
#define problemDIM 2
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
//n, c
#define numFields 2

//define time step parameters
#define timeStep 1.0e-4
#define timeFinal 1.0
#define timeIncrements 100000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 10
#define outputFileDirectory "."

//define Cahn-Hilliard parameters (No Gradient energy term)
#define McV 1.0

//define Allen-Cahn parameters
#define MnV 1.0
#define KnV 1.0

//define free energy expressions
#define faV (-1.6704-4.776*c+5.1622*c*c-2.7375*c*c*c+1.3687*c*c*c*c)
#define facV (-4.776 + 10.3244*c - 8.2125*c*c + 5.4748*c*c*c)
#define faccV (10.3244-16.425*c+16.4244*c*c)
#define fbV (5.0*c*c-5.9746*c-1.5924)
#define fbcV (10.0*c-5.9746)
#define fbccV (10.0)
#define hV (10.0*n*n*n-15.0*n*n*n*n+6.0*n*n*n*n*n)
#define hnV (30.0*n*n-60.0*n*n*n+30.0*n*n*n*n)

//define required residuals
//#define rcxV ( cx*((1.0-hV)*faccV+hV*fbccV) + nx*((fbcV-facV)*hnV) )
//#define rnV  ( (fbV-faV)*hnV )
//#define rnxV nx
#define muxV ( cx*((1.0-hV)*faccV+hV*fbccV) + nx*((fbcV-facV)*hnV) )
#define rcV   (c)
#define rcxV  (constV(-McV*timeStep)*muxV)
#define rnV  (n-constV(timeStep*MnV)*(fbV-faV)*hnV)
#define rnxV (constV(-timeStep*KnV*MnV)*nx)

