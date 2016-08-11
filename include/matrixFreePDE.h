//base class for matrix Free implementation of PDE's
#ifndef MATRIXFREEPDE_H
#define MATRIXFREEPDE_H

//general headers
#include <fstream>
#include <sstream>
#include <iterator> // is this necessary?

//dealii headers
#include "dealIIheaders.h"
#include "defaultValues.h"

//PRISMS headers
#include "fields.h"


//define data types
#ifndef scalarType
typedef dealii::VectorizedArray<double> scalarType;
#endif
#ifndef vectorType
typedef dealii::parallel::distributed::Vector<double> vectorType;
#endif
//define FE system types
#ifndef typeScalar
typedef dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,1,double>           typeScalar;
#endif
#ifndef typeVector
typedef dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,problemDIM,double>  typeVector;
#endif
//define data value types
#ifndef scalarvalueType
typedef dealii::VectorizedArray<double> scalarvalueType;
#endif
#ifndef vectorvalueType
typedef dealii::Tensor<1, problemDIM, dealii::VectorizedArray<double> > vectorvalueType;
#endif
#if problemDIM==1
#ifndef scalargradType
typedef dealii::VectorizedArray<double> scalargradType;
#endif
#ifndef vectorgradType
typedef dealii::VectorizedArray<double> vectorgradType;
#endif
#ifndef vectorhessType
typedef dealii::VectorizedArray<double> vectorhessType;
#endif
#else
#ifndef scalargradType
typedef dealii::Tensor<1, problemDIM, dealii::VectorizedArray<double> > scalargradType;
#endif
#ifndef scalarhessType
typedef dealii::Tensor<2,problemDIM,dealii::VectorizedArray<double> > scalarhessType;
#endif
#ifndef vectorgradType
typedef dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > vectorgradType;
#endif
#ifndef vectorhessType
typedef dealii::Tensor<3, problemDIM, dealii::VectorizedArray<double> > vectorhessType;
#endif
#endif

#include "model_variables.h"

//macro for constants
#define constV(a) make_vectorized_array(a)
//macro for defining subdomain specific functions
#define subdomain(geometricExpression, functionExpression)  ( (geometricExpression) ? (functionExpression) : constV(0.0))

//
using namespace dealii;
//
//base class for matrix free PDE's
//
/**
 * This is the base class for matrix free implementation for Parabolic and Elliptic BVP's.
 * This class supports MPI, Threads and Vectorization (Hybrid Parallel). 
 */
template <int dim>
class MatrixFreePDE:public Subscriptor
{
 public:
  /**
   * Class contructor
   */
  MatrixFreePDE(); 
  ~MatrixFreePDE(); 

  /**
   * Initializes the data structures
   */
  void init  (unsigned int iter=0);
  void initForTests();

  /**
   * Solve's the system of equations
   */
  void solve ();
  void vmult (vectorType &dst, const vectorType &src) const;
  std::vector<Field<dim> >                  fields;
  //matrix free methods
  void computeInvM();
  
  // Virtual function to shift the concentration
  void shiftConcentration();

 protected:
  void solveIncrement ();
  void outputResults  ();
  parallel::distributed::Triangulation<dim> triangulation;
  std::vector<FESystem<dim>*>          FESet;
  std::vector<const ConstraintMatrix*> constraintsSet, constraintsHangingNodesSet;
  std::vector<const DoFHandler<dim>*>  dofHandlersSet;
  std::vector<const IndexSet*>         locally_relevant_dofsSet;
  std::vector<ConstraintMatrix*>       constraintsSet2, constraintsHangingNodesSet2;
  std::vector<DoFHandler<dim>*>        dofHandlersSet2;
  std::vector<IndexSet*>               locally_relevant_dofsSet2;
  std::vector<vectorType*>             solutionSet, residualSet;
  std::vector<parallel::distributed::SolutionTransfer<dim, vectorType>*> soltransSet;
  
  //matrix free objects
  MatrixFree<dim,double>               matrixFreeObject;
  vectorType                           invM, dU;
  
  //matrix free methods
  unsigned int currentFieldIndex;
  unsigned int num_quadrature_points;
  void computeRHS();

  //AMR methods
  void refineGrid();
  void refineMesh(unsigned int _currentIncrement);
  virtual void adaptiveRefine(unsigned int _currentIncrement);
  virtual void adaptiveRefineCriterion();
  
  //virtual methods to be implemented in the derived class
  //method to calculate LHS(implicit)
  virtual void getLHS(const MatrixFree<dim,double> &data, 
		      vectorType &dst, 
		      const vectorType &src,
		      const std::pair<unsigned int,unsigned int> &cell_range) const;
  //method to calculate RHS (implicit/explicit) -- abstract method
  virtual void getRHS (const MatrixFree<dim,double> &data, 
		       std::vector<vectorType*> &dst, 
		       const std::vector<vectorType*> &src,
		       const std::pair<unsigned int,unsigned int> &cell_range) const = 0;
  
  //methods to apply dirichlet BC's
  std::vector<std::map<dealii::types::global_dof_index, double>*> valuesDirichletSet;
  virtual void markBoundaries();
  virtual void applyDirichletBCs();

  //methods to apply initial conditions
  virtual void applyInitialConditions();
  virtual void modifySolutionFields ();
 
  void computeEnergy();
  virtual void getEnergy(const MatrixFree<dim,double> &data,
		    std::vector<vectorType*> &dst,
		    const std::vector<vectorType*> &src,
		    const std::pair<unsigned int,unsigned int> &cell_range);

  //utility functions
  //return index of given field name if exists, else throw error
  unsigned int getFieldIndex(std::string _name);


  std::vector<double> freeEnergyValues;
  void outputFreeEnergy(std::vector<double>& freeEnergyValues);

  // Virtual method to compute the integral of a field
  void computeIntegral(double& integratedField);

  //variables for time dependent problems 
  //isTimeDependentBVP flag is used to see if invM, time steppping in
  //run(), etc are necessary
  bool isTimeDependentBVP, isEllipticBVP;
  unsigned int parabolicFieldIndex, ellipticFieldIndex;
  double dtValue, currentTime, finalTime;
  unsigned int currentIncrement, totalIncrements;

  //parallel message stream
  ConditionalOStream  pcout;  
  //compute time log
  mutable TimerOutput computing_timer;

  double energy;
  std::vector<double> energy_components;
};

//other matrixFree headers 
//(these are source files, which will are temporarily treated as
//header files till library packaging scheme is finalized)
#include "../src/matrixfree/matrixFreePDE.cc"
#include "../src/matrixfree/init.cc"
#include "../src/matrixfree/initForTests.cc"
#include "../src/matrixfree/refine.cc"
#include "../src/matrixfree/invM.cc"
#include "../src/matrixfree/computeLHS.cc"
#include "../src/matrixfree/computeRHS.cc"
#include "../src/matrixfree/modifyFields.cc"
#include "../src/matrixfree/solve.cc"
#include "../src/matrixfree/solveIncrement.cc"
#include "../src/matrixfree/outputResults.cc"
#include "../src/matrixfree/markBoundaries.cc"
#include "../src/matrixfree/boundaryConditions.cc"
#include "../src/matrixfree/initialConditions.cc"
#include "../src/matrixfree/utilities.cc"
#include "../src/matrixfree/calcFreeEnergy.cc"
#include "../src/matrixfree/integrate_and_shift_field.cc"

#endif
