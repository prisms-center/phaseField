//base class for matrix Free implementation of PDE's
#ifndef MATRIXFREEPDE_H
#define MATRIXFREEPDE_H

//general headers
#include <fstream>
#include <sstream>

//dealii headers
#include "dealIIheaders.h"
#include "defaultValues.h"

//PRISMS headers
#include "fields.h"

//define data types  
typedef dealii::parallel::distributed::Vector<double> vectorType;
typedef dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,1,double>           typeScalar;
typedef dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,problemDIM,double>  typeVector;
//
#if problemDIM==1
typedef dealii::VectorizedArray<double> valueType;
typedef dealii::VectorizedArray<double> scalarvalueType;
typedef dealii::VectorizedArray<double> gradType;
typedef dealii::VectorizedArray<double> scalargradType;
#else 
typedef dealii::Tensor<1, numFields, dealii::VectorizedArray<double> > valueType;
typedef dealii::VectorizedArray<double> scalarvalueType;
typedef dealii::Tensor<1, numFields, dealii::Tensor<1, problemDIM, dealii::VectorizedArray<double> > > gradType;
typedef dealii::Tensor<1, problemDIM, dealii::VectorizedArray<double> > scalargradType;
#endif
typedef dealii::VectorizedArray<double> scalarType;

//macro for constants
#define constV(a) make_vectorized_array(a)

//
using namespace dealii;
//
//base class for matrix free PDE's
//
/**
 * This class collects all the data that is stored for the matrix free
 * implementation. The storage scheme is tailored towards several loops
 * performed with the same data, i.e., typically doing many matrix-vector
 * products or residual computations on the same mesh.
 */
template <int dim>
class MatrixFreePDE:public Subscriptor
{
 public:
  /**
   * Collects the options for initialization of the MatrixFree class. The
   * first parameter specifies the MPI communicator to be used, the second the
   * parallelization options in shared memory (task-based parallelism, where
   * one can choose between no parallelism and three schemes that avoid that
   * cells with access to the same vector entries are accessed
   * simultaneously), the third with the block size for task parallel
   * scheduling, the fourth the update flags that should be stored by this
   * class.
   *
   * The fifth parameter specifies the level in the triangulation from which
   * the indices are to be used. If the level is set to
   * numbers::invalid_unsigned_int, the active cells are traversed, and
   * otherwise the cells in the given level. This option has no effect in case
   * a DoFHandler or hp::DoFHandler is given.
   *
   * The parameter @p initialize_plain_indices indicates whether the DoFInfo
   * class should also allow for access to vectors without resolving
   * constraints.
   *
   * The last two parameters allow the user to disable some of the
   * initialization processes. For example, if only the scheduling that avoids
   * touching the same vector/matrix indices simultaneously is to be found,
   * the mapping needs not be initialized. Likewise, if the mapping has
   * changed from one iteration to the next but the topology has not (like
   * when using a deforming mesh with MappingQEulerian), it suffices to
   * initialize the mapping only.
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
  vectorType                           invM;
  
  //matrix free methods
  unsigned int currentFieldIndex;
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

#endif
