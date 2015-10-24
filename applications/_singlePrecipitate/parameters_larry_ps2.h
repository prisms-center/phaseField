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
#define dt 1.67e-5
#define numIncrements 7000000

//define solver tolerances
#define relativeResTol 1.0e-10
#define absoluteResTol 1.0e-14
#define dUTol 1.0e-10

//define results output parameters
#define writeOutput true
#define skipOutputSteps 1000000

//define general physics parameters
#define numStructuralOrderParameters 3

//define Cahn-Hilliard parameters (No Gradient energy)
#define Mc 1.0

//define Allen-Cahn parameters
#define MnVals {50.0,50.0,50.0}
//#define MnVals {50.0}
#define Kn0Tensor {0.0150,0,0,0,0.0188,0,0,0,0.00571}
#define Kn1Tensor {0.123,0,0,0,0.123,0,0,0,0.123}
#define Kn2Tensor {0.123,0,0,0,0.123,0,0,0,0.123}
//define energy barrier coefficient
#define W -1.0

//define Mechanical properties
#define MaterialModelv ISOTROPIC
#define MaterialConstantsv {1.0,0.3}
#define sf0StrainV {0,0,0,0,0,0,0,0,0}
#define sf1StrainV {0,0,0,0,0,0,0,0,0}
#define sf2StrainV {0,0,0,0,0,0,0,0,0}
#define skipElasticitySteps 100000000

//define free energy expressions
#define faV (24.7939*c*c - 1.6752*c - 1.9453e-06)
#define facV (49.5878*c - 1.6752)
#define faccV (49.5878)
#define fbV (37.9316*c*c - 10.7373*c + 0.5401)
#define fbcV (75.8633*c - 10.7373)
#define fbccV (75.8633)
#define h0V (3.0*n[0]*n[0]-2.0*n[0]*n[0]*n[0])
#define h1V (3.0*n[1]*n[1]-2.0*n[1]*n[1]*n[1])
#define h2V (3.0*n[2]*n[2]-2.0*n[2]*n[2]*n[2])
#define h0nV (6.0*n[0]-6.0*n[0]*n[0])
#define h1nV (6.0*n[1]-6.0*n[1]*n[1])
#define h2nV (6.0*n[2]-6.0*n[2]*n[2])
#define fbarrierV (n[0]*n[0]-2.0*n[0]*n[0]*n[0]+n[0]*n[0]*n[0]*n[0])
#define fbarriernV (2.0*n[0]-6.0*n[0]*n[0]+4.0*n[0]*n[0]*n[0])

//define required residuals
#define rcxV ( cx*((1.0-h0V-h1V-h2V)*faccV+(h0V+h1V+h2V)*fbccV) + nx[0]*((fbcV-facV)*h0nV) + nx[1]*((fbcV-facV)*h1nV) + nx[2]*((fbcV-facV)*h2nV) )
#define rn0V  ( (fbV-faV)*h0nV + W*fbarriernV)
#define rn1V  ( (fbV-faV)*h1nV )
#define rn2V  ( (fbV-faV)*h2nV )
#define rn0xV nx[0]
#define rn1xV nx[1]
#define rn2xV nx[2]
