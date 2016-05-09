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

