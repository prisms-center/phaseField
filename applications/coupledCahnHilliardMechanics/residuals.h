//define number of fields in the problem
//c, mu, u
#define numFields (2+problemDIM)

//define Cahn-Hilliard parameters
#define McV 1.0
#define KcV 1.0

//define material properties
#define MaterialModelV ISOTROPIC
#define MaterialConstantsV {1.0,0.3}

//define free energy expressions
#define fcV (4.0*c*(c-1.0)*(c-0.5))

//define required residuals
#define rmuV  (fcV)
#define rmuxV (constV(KcV)*cx)
#define rcV   (c)
#define rcxV  (constV(-McV*timeStep)*mux)

//define chemical strain
#define chemicalStrainV (-0.01*c)
