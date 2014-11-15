//Parameters list for beta prime precipitation evolution problem

//define problem dimensions
#define problemDIM 2
#define spanX 100.0
#define spanY 100.0
#define spanZ 100.0

//define mesh parameters
#define refineFactor 8
#define finiteElementDegree 1

//define time step parameters
#define dt 1.0e-4
#define numIncrements 20000

//define solver tolerances
#define relativeResTol 1.0e-10
#define absoluteResTol 1.0e-14
#define dUTol 1.0e-10

//define results output parameters
#define writeOutput true
#define skipOutputSteps 2000

//define general physics parameters
#define numStructuralOrderParameters 3

//define Cahn-Hilliard parameters (No Gradient energy)
#define Mc 1.0

//define Allen-Cahn parameters
#define MnVals {1.0,1.0,1.0}
#define Kn0Tensor {{4.0,0,0},{0,1.0,0},{0,0,1.0}}
#define Kn1Tensor {{1.75,-1.299,0},{-1.299,3.25,0},{0,0,1.0}}
#define Kn2Tensor {{1.75, 1.299,0},{1.299,3.25,0},{0,0,1.0}}

//define Mechanical properties
#define MaterialModelv ISOTROPIC
#define MaterialConstantsv {1.0,0.3}
#define sf0StrainV {{0.0345,0,0},{0,0.0185,0},{0,0,-0.00270}}
#define sf1StrainV {{0.0225,-0.0069,0},{-0.0069,0.0305,0},{0,0,-0.00270}}
#define sf2StrainV {{0.0225, 0.0069,0},{0.0069,0.0305,0},{0,0,-0.00270}}
#define skipElasticitySteps 100

//define free energy expressions
#define faV (-1.6704-4.776*c+5.1622*c*c-2.7375*c*c*c+1.3687*c*c*c*c)
#define facV (-4.776 + 10.3244*c - 8.2125*c*c + 5.4748*c*c*c)
#define faccV (10.3244-16.425*c+16.4244*c*c)
#define fbV (5.0*c*c-5.9746*c-1.5924)
#define fbcV (10.0*c-5.9746)
#define fbccV (10.0)
#define h0V (10.0*n[0]*n[0]*n[0]-15.0*n[0]*n[0]*n[0]*n[0]+6.0*n[0]*n[0]*n[0]*n[0]*n[0])
#define h1V (10.0*n[1]*n[1]*n[1]-15.0*n[1]*n[1]*n[1]*n[1]+6.0*n[1]*n[1]*n[1]*n[1]*n[1])
#define h2V (10.0*n[2]*n[2]*n[2]-15.0*n[2]*n[2]*n[2]*n[2]+6.0*n[2]*n[2]*n[2]*n[2]*n[2])
#define h0nV (30.0*n[0]*n[0]-60.0*n[0]*n[0]*n[0]+30.0*n[0]*n[0]*n[0]*n[0])
#define h1nV (30.0*n[1]*n[1]-60.0*n[1]*n[1]*n[1]+30.0*n[1]*n[1]*n[1]*n[1])
#define h2nV (30.0*n[2]*n[2]-60.0*n[2]*n[2]*n[2]+30.0*n[2]*n[2]*n[2]*n[2])

//define required residuals
#define rcxV ( cx*((1.0-h0V-h1V-h2V)*faccV+(h0V+h1V+h2V)*fbccV) + nx[0]*((fbcV-facV)*h0nV) + nx[1]*((fbcV-facV)*h1nV) + nx[2]*((fbcV-facV)*h2nV) )
#define rn0V  ( (fbV-faV)*h0nV )
#define rn1V  ( (fbV-faV)*h1nV )
#define rn2V  ( (fbV-faV)*h2nV )
#define rn0xV nx[0]
#define rn1xV nx[1]
#define rn2xV nx[2]

