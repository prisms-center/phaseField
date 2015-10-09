//Parameters list for beta prime precipitation evolution problem

//define problem dimensions
#define problemDIM 3
#define spanX 16.0
#define spanY 16.0
#define spanZ 8.0

//define mesh parameters
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define refineFactor 7
#define finiteElementDegree 1

//define time step parameters
#define timeStep 1.67e-5
#define timeFinal 1588836  // This should really be tied to an expression, not directly defined
#define timeIncrements 1000000 //7000000
#define skipImplicitSolves 1000

//define solver paramters
#define solverType SolverCG
#define relSolverTolerance 1.0e-10
#define maxSolverIterations 1000

//define results output parameters
#define writeOutput true
#define skipOutputSteps 100000

#define numFields (4+problemDIM)

//define Cahn-Hilliard parameters (No Gradient energy)
#define McV 1.0

//define Allen-Cahn parameters
#define Mn1V 1.0
#define Mn2V 1.0
#define Mn3V 1.0

double Kn1[3][3]={{0.0150,0,0},{0,0.0188,0},{0,0,0.00571}};
double Kn2[3][3]={{0.123,0,0},{0,0.123,0},{0,0,0.123}};
double Kn3[3][3]={{0.123,0,0},{0,0.123,0},{0,0,0.123}};

//define energy barrier coefficient
#define W -1.0


//define Mechanical properties
#define MaterialModelV ISOTROPIC
#define MaterialConstantsV {1.0,0.3}
double sf1Strain[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sf2Strain[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sf3Strain[3][3] = {{0,0,0},{0,0,0},{0,0,0}};


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

/* Parameters from Larry's code
//define required residuals
#define rcxV ( cx*((1.0-h0V-h1V-h2V)*faccV+(h0V+h1V+h2V)*fbccV) + nx[0]*((fbcV-facV)*h0nV) + nx[1]*((fbcV-facV)*h1nV) + nx[2]*((fbcV-facV)*h2nV) )
#define rn1V  ( (fbV-faV)*h0nV + W*fbarriernV)
#define rn2V  ( (fbV-faV)*h1nV )
#define rn3V  ( (fbV-faV)*h2nV )
#define rn1xV nx[0]
#define rn2xV nx[1]
#define rn3xV nx[2]
*/
#define rcV   (c)
#define rcxTemp ( cx*((1.0-h1V-h2V-h3V)*faccV+(h1V+h2V+h3V)*fbccV) + n1x*((fbcV-facV)*hn1V) + n2x*((fbcV-facV)*hn2V) + n3x*((fbcV-facV)*hn3V) )
#define rcxV  (constV(-timeStep*McV)*rcxTemp)

#define rn1V   (n1-constV(-timeStep*Mn1V)*((fbV-faV)*hn1V+W*fbarriernV-CEE1))
#define rn2V   (n2-constV(-timeStep*Mn2V)*((fbV-faV)*hn2V+W*fbarriernV-CEE2))
#define rn3V   (n3-constV(-timeStep*Mn3V)*((fbV-faV)*hn3V+W*fbarriernV-CEE3))
#define rn1xV  (constV(-timeStep*Mn1V)*Knx1)
#define rn2xV  (constV(-timeStep*Mn2V)*Knx2)
#define rn3xV  (constV(-timeStep*Mn3V)*Knx3)

/* Parameters from the CH_AC_mech example code
//define Mechanical properties
#define MaterialModelV ISOTROPIC
#define MaterialConstantsV {1.0,0.3}
double sf1Strain[3][3]={{0.0345,0,0},{0,0.0185,0},{0,0,-0.00270}};
double sf2Strain[3][3]={{0.0225,-0.0069,0},{-0.0069,0.0305,0},{0,0,-0.00270}};
double sf3Strain[3][3]={{0.0225, 0.0069,0},{0.0069,0.0305,0},{0,0,-0.00270}};


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


//define required residuals
#define rcV   (c)
#define rcxTemp ( cx*((1.0-h1V-h2V-h3V)*faccV+(h1V+h2V+h3V)*fbccV) + n1x*((fbcV-facV)*hn1V) + n2x*((fbcV-facV)*hn2V) + n3x*((fbcV-facV)*hn3V) )
#define rcxV  (constV(-timeStep*McV)*rcxTemp)
//n1
#define rn1V   (n1-constV(-timeStep*Mn1V)*((fbV-faV)*hn1V-CEE1))
#define rn2V   (n2-constV(-timeStep*Mn2V)*((fbV-faV)*hn2V-CEE2))
#define rn3V   (n3-constV(-timeStep*Mn3V)*((fbV-faV)*hn3V-CEE3))
#define rn1xV  (constV(-timeStep*Mn1V)*Knx1)
#define rn2xV  (constV(-timeStep*Mn2V)*Knx2)
#define rn3xV  (constV(-timeStep*Mn3V)*Knx3)

*/



