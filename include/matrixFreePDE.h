//base class for matrix Free implementation of PDE's
#ifndef MATRIXFREEPDE_H
#define MATRIXFREEPDE_H

//general headers
#include <fstream>
#include <sstream>

//dealii headers
#include "dealIIheaders.h"

//PRISMS headers
#include "fields.h"

//define data types  
typedef dealii::parallel::distributed::Vector<double> vectorType;
typedef dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,1,double>           typeScalar;
typedef dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,problemDIM,double>  typeVector;
#if problemDIM==1
typedef dealii::VectorizedArray<double> gradType;
#else 
typedef dealii::Tensor<1, problemDIM, dealii::VectorizedArray<double> > gradType;
#endif

using namespace dealii;

//
//base class for matrix free PDE's
//
template <int dim>
class MatrixFreePDE:public Subscriptor
{
 public:
  MatrixFreePDE(); 
  ~MatrixFreePDE(); 
  void init  ();
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
  vectorType                           invM;
  
  //matrix free methods
  void computeInvM();
  void updateRHS();
  void computeRHS(const MatrixFree<dim,double> &data, 
		  std::vector<vectorType*> &dst, 
		  const std::vector<vectorType*> &src,
		  const std::pair<unsigned int,unsigned int> &cell_range) const;
  void computeLHS(const MatrixFree<dim,double> &data, 
		  vectorType &dst, 
		  const vectorType &src,
		  const std::pair<unsigned int,unsigned int> &cell_range) const;

  //virtual method to be implemented in derived classe
  virtual void getRHS(std::map<std::string, typeScalar*>  valsScalar,	\
		      std::map<std::string, typeVector*>  valsVector,	\
		      unsigned int q) const = 0;  
  
  //residual, matrix-vector computation variables
  std::map<std::string,bool> getValue, setValue;
  std::map<std::string,bool> getGradient, setGradient;

  //variables for time dependent problems 
  //isTimeDependentBVP flag is used to see if invM, time steppping in
  //run(), etc are necessary
  bool isTimeDependentBVP;
  double timeStep, currentTime, finalTime;
  unsigned int currentIncrement, totalIncrements;
  
  //parallel message stream
  ConditionalOStream  pcout;  
  //compute time log
  TimerOutput computing_timer;
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


#endif
