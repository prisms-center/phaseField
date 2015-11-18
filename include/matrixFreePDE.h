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
typedef dealii::VectorizedArray<double> scalarType;
typedef dealii::parallel::distributed::Vector<double> vectorType;
//define FE system types
typedef dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,1,double>           typeScalar;
typedef dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,problemDIM,double>  typeVector;
//define data value types
typedef dealii::VectorizedArray<double> scalarvalueType;
#if problemDIM==1
typedef dealii::VectorizedArray<double> scalargradType;
typedef dealii::VectorizedArray<double> vectorgradType;
#else 
typedef dealii::Tensor<1, problemDIM, dealii::VectorizedArray<double> > scalargradType;
typedef dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > vectorgradType;
#endif

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
  void init  ();

  /**
   * Solve's the system of equations
   */
  void solve ();
  void vmult (vectorType &dst, const vectorType &src) const;
  std::vector<Field<dim> >                  fields;

 protected:
  void solveIncrement ();
  void outputResults  ();
  parallel::distributed::Triangulation<dim> triangulation;
  std::vector<FESystem<dim>*>          FESet;
  std::vector<const ConstraintMatrix*> constraintsSet;
  std::vector<const DoFHandler<dim>*>  dofHandlersSet;
  std::vector<const IndexSet*>         locally_relevant_dofsSet;
  std::vector<vectorType*>             solutionSet, residualSet;

  //matrix free objects
  MatrixFree<dim,double>               matrixFreeObject;
  vectorType                           invM, dU;
  
  //matrix free methods
  unsigned int currentFieldIndex;
  unsigned int num_quadrature_points;
  void computeInvM();
  //RHS
  void updateRHS();
  void computeRHS(const MatrixFree<dim,double> &data, 
		  std::vector<vectorType*> &dst, 
		  const std::vector<vectorType*> &src,
		  const std::pair<unsigned int,unsigned int> &cell_range) const;
  //virtual methods to be implemented in derived classe
  //methods to calculate RHS (implicit/explicit)
  virtual void getRHS(std::map<std::string, typeScalar*>  valsScalar,	\
		      std::map<std::string, typeVector*>  valsVector,	\
		      unsigned int q) const = 0;  
  //LHS
  template <typename T>
    void computeLHS(const MatrixFree<dim,double> &data, 
		    vectorType &dst, 
		    const vectorType &src,
		    const std::pair<unsigned int,unsigned int> &cell_range) const;
  //methods to calculate LHS(implicit)
  virtual void getLHS(typeScalar& vals, unsigned int q) const;  
  virtual void getLHS(typeVector& vals, unsigned int q) const;  
  
  //methods to apply dirichlet BC's
  virtual void markBoundaries();
  virtual void applyDirichletBCs();
  virtual void applyInitialConditions();

  //utility functions
  //return index of given field name if exists, else throw error
  unsigned int getFieldIndex(std::string _name);
  virtual void computeFreeEnergyValue(std::vector<double>& freeEnergyValues);
  std::vector<double> freeEnergyValues;
  void outputFreeEnergy(std::vector<double>& freeEnergyValues);

  //residual, matrix-vector computation variables
  std::map<std::string,bool> getValue, setValue;
  std::map<std::string,bool> getGradient, setGradient;

  //variables for time dependent problems 
  //isTimeDependentBVP flag is used to see if invM, time steppping in
  //run(), etc are necessary
  bool isTimeDependentBVP;
  double dtValue, currentTime, finalTime;
  unsigned int currentIncrement, totalIncrements;

  //parallel message stream
  ConditionalOStream  pcout;  
  //compute time log
  mutable TimerOutput computing_timer;
};

//other matrixFree headers 
//(these are source files, which will are temporarily treated as
//header files till library packaging scheme is finalized)
#include "../src/matrixfree/matrixFreePDE.cc"
#include "../src/matrixfree/init.cc"
#include "../src/matrixfree/invM.cc"
#include "../src/matrixfree/computeLHS.cc"
#include "../src/matrixfree/computeRHS.cc"
#include "../src/matrixfree/solve.cc"
#include "../src/matrixfree/solveIncrement.cc"
#include "../src/matrixfree/outputResults.cc"
#include "../src/matrixfree/markBoundaries.cc"
#include "../src/matrixfree/boundaryConditions.cc"
#include "../src/matrixfree/initialConditions.cc"
#include "../src/matrixfree/utilities.cc"
#include "../src/matrixfree/freeEnergy.cc"

#endif
