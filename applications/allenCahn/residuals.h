//define Allen-Cahn parameters
#define MnV 1.0
#define KnV 4.0

//define the derivative of the free energy with respect to n
#define fnV (4.0*n*(n-1.0)*(n-0.5))

//define required residuals
#define rnV  (n-constV(timeStep*MnV)*fnV)
#define rnxV (constV(-timeStep*KnV*MnV)*nx)
