//some default values
#ifndef DEFAULTS_H
#define DEFAULTS_H

//problem dimension (default value:3)
#ifndef problemDIM
#define problemDIM 3
#endif

//spanX (default value:1.0)
#ifndef spanX
#define spanX 1.0
#endif

//spanX (default value:1.0)
#ifndef spanY
#define spanY 1.0
#endif

//spanX (default value:1.0)
#ifndef spanZ
#define spanZ 1.0
#endif

//FE mesh refine factor (default value:3, which gives 8^problemDIM elements)
#ifndef refineFactor
#define refineFactor 3
#endif

//FE polynomial order (default value:1)
#ifndef finiteElementDegree
#define finiteElementDegree 1
#endif

//number of fields (default value:1)
#ifndef numFields
#define numFields 1
#endif

//write output files (default value:true)
#ifndef writeOutput 
#define writeOutput true
#endif

//output written at every n'th increment for time dependent problems (default value:1)
#ifndef skipOutputSteps
#define skipOutputSteps 1
#endif

//solver type for implcit solves (default value:SolverCG)
#ifndef solverType
#define solverType SolverCG
#endif

//relative linear solver tolerance for implicit solves (default value:1.0e-10)
#ifndef solverTolerance
#define solverTolerance 1.0e-10
#endif

//max number of linear solver iterations for implicit solves (default value:1000)
#ifndef maxSolverIterations
#define maxSolverIterations 1000
#endif

//number of implicit solves to skip. None are skipped if value is 1, which is the default.
#ifndef skipImplicitSolves
#define skipImplicitSolves 1
#endif

//number of implicit solves to skip. None are skipped if value is 1, which is the default.
#ifndef outputCondition
#define outputCondition "EQUAL_SPACING"
#endif

#ifndef outputList
#define outputList {0}
#endif

#endif
