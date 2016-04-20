//Parameters list for the coupled Cahn Hilliard and Allen Cahn problem

//define problem dimensions
#define problemDIM 3
#define spanX 64.0
#define spanY 64.0
#define spanZ 64.0

//define mesh parameters
#define subdivisionsX 5 //1 //5
#define subdivisionsY 5 //1 //5
#define subdivisionsZ 5 //1 //5
#define refineFactor 3 //7 //3
#define finiteElementDegree 2 //1 //2

//define number of fields in the problem
//n, c
#define numFields 2

//define time step parameters
#define timeStep (5.0/5500.0) //5.0e-4 //(5.0/5500.0) //5.0e-4
#define timeFinal 10000.0
#define timeIncrements 27500 //50000 //27500 //50000 

//define results output parameters
#define writeOutput true
#define skipOutputSteps 50000

//define Cahn-Hilliard parameters (No Gradient energy term)
#define McV 1.0

//define Allen-Cahn parameters
#define MnV 5.0
#define KnV 0.5

//define free energy expressions
#define faV (24.7939*c*c - 1.6752*c - 1.9453e-06)
#define facV (49.5878*c - 1.6752)
#define faccV (49.5878)
#define fbV (37.9316*c*c - 10.7373*c + 0.5401)
#define fbcV (75.8633*c - 10.7373)
#define fbccV (75.8633)
#define hV (3.0*n*n-2.0*n*n*n)
#define hnV (6.0*n-6.0*n*n)

//define required residuals
//#define rcxV ( cx*((1.0-hV)*faccV+hV*fbccV) + nx*((fbcV-facV)*hnV) )
//#define rnV  ( (fbV-faV)*hnV )
//#define rnxV nx
#define muxV ( cx*((1.0-hV)*faccV+hV*fbccV) + nx*((fbcV-facV)*hnV) )
#define rcV   (c)
#define rcxV  (constV(-McV*timeStep)*muxV)
#define rnV  (n-constV(timeStep*MnV)*(fbV-faV)*hnV)
#define rnxV (constV(-timeStep*KnV*MnV)*nx)

#define avg_c 0.02
#define initial_radius 1.0
#define x_denom 36.0
#define y_denom 36.0
#define z_denom 16.0
#define initial_interface_coeff 0.5

