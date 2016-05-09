// Define Cahn-Hilliard parameters (no gradient energy terms)
#define McV 1.0

// Define Allen-Cahn parameters
#define Mn1V 1.0
#define Mn2V 1.0
#define Mn3V 1.0
double Kn1[3][3]={{4.0,0,0},{0,1.0,0},{0,0,1.0}};
double Kn2[3][3]={{1.75,-1.299,0},{-1.299,3.25,0},{0,0,1.0}};
double Kn3[3][3]={{1.75, 1.299,0},{1.299,3.25,0},{0,0,1.0}};

// Define Mechanical properties
#define n_dependent_stiffness false
// Mechanical symmetry of the material and stiffness parameters
// Used throughout system if n_dependent_stiffness == false, used in n=0 phase if n_dependent_stiffness == true
#define MaterialModelV ISOTROPIC
#define MaterialConstantsV {2.0,0.3}

// Used in n=1 phase if n_dependent_stiffness == true
#define MaterialModelBetaV ISOTROPIC
#define MaterialConstantsBetaV {2.0,0.3}

// Stress-free transformation strains
// Linear fits for the stress-free transformation strains in for sfts = sfts_linear * c + sfts_const
double sfts_linear1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sfts_const1[3][3] = {{0.0345,0,0},{0,0.0185,0},{0,0,-0.00270}};

double sfts_linear2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sfts_const2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

double sfts_linear3[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sfts_const3[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

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

// Define required residuals
#define rcV   (c)
#define rcxTemp ( cx*((1.0-h1V-h2V-h3V)*faccV+(h1V+h2V+h3V)*fbccV) + n1x*((fbcV-facV)*hn1V) + n2x*((fbcV-facV)*hn2V) + n3x*((fbcV-facV)*hn3V) + grad_mu_el)
#define rcxV  (constV(-timeStep*McV)*rcxTemp)

#define rn1V   (n1-constV(timeStep*Mn1V)*((fbV-faV)*hn1V+nDependentMisfitAC1+heterMechAC1))
#define rn2V   (n2-constV(timeStep*Mn2V)*((fbV-faV)*hn2V+nDependentMisfitAC2+heterMechAC2))
#define rn3V   (n3-constV(timeStep*Mn3V)*((fbV-faV)*hn3V+nDependentMisfitAC3+heterMechAC3))
#define rn1xV  (constV(-timeStep*Mn1V)*Knx1)
#define rn2xV  (constV(-timeStep*Mn2V)*Knx2)
#define rn3xV  (constV(-timeStep*Mn3V)*Knx3)
