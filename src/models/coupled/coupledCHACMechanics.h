//Matrix Free implementation of coupled Cahn-Hilliard, Allen-Cahn and Mechanics formulation 
#ifndef CHACMECHANICS_H
#define CHACMECHANICS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

//material models
#include "../mechanics/computeStress.h"

template <int dim>
class CoupledCHACMechanicsProblem: public MatrixFreePDE<dim>
{
 public: 
  CoupledCHACMechanicsProblem();

 private:
  //elasticity matrix
  Table<2, double> CIJ;

  //RHS implementation for explicit solve
  void getRHS(std::map<std::string, typeScalar*>  valsScalar, \
	      std::map<std::string, typeVector*>  valsVector, \
	      unsigned int q) const;

  //LHS implementation for implicit solve 
  void getLHS(typeVector& vals, unsigned int q) const;  

  //method to apply initial conditions
  void applyInitialConditions();
 
  //methods to apply dirichlet BC's on displacement
  void applyDirichletBCs();


  // method to modify the fields for nucleation
  void modifySolutionFields();

  // calculate (and output) the total free energy of the system
  void computeFreeEnergyValue(std::vector<double>& freeEnergyValues);

};

//constructor
template <int dim>
CoupledCHACMechanicsProblem<dim>::CoupledCHACMechanicsProblem(): MatrixFreePDE<dim>(),
  CIJ(2*dim-1+dim/3,2*dim-1+dim/3)
{
  //initialize elasticity matrix
#if defined(MaterialModelV) && defined(MaterialConstantsV)
  double materialConstants[]=MaterialConstantsV;
  getCIJMatrix<dim>(MaterialModelV, materialConstants, CIJ, this->pcout);
#else
#error Compile ERROR: missing material property variable: MaterialModelV, MaterialConstantsV
#endif

  //"c"
  this->getValue["c"]=true; this->getGradient["c"]=true;
  this->setValue["c"]=true; this->setGradient["c"]=true;
  //"n1"
  this->getValue["n1"]=true; this->getGradient["n1"]=true;
  this->setValue["n1"]=true; this->setGradient["n1"]=true;
 //"n2"
  this->getValue["n2"]=true; this->getGradient["n2"]=true;
  this->setValue["n2"]=true; this->setGradient["n2"]=true;
 //"n3"
  this->getValue["n3"]=true; this->getGradient["n3"]=true;
  this->setValue["n3"]=true; this->setGradient["n3"]=true;
  //"u"
  this->getValue["u"]=false; this->getGradient["u"]=true;
  this->setValue["u"]=false; this->setGradient["u"]=true;
}


template <int dim>
void  CoupledCHACMechanicsProblem<dim>::getRHS(std::map<std::string, typeScalar*>  valsScalar, \
					       std::map<std::string, typeVector*>  valsVector, \
					       unsigned int q) const{
  //"c" fields
  scalarvalueType c = valsScalar["c"]->get_value(q);
  scalargradType cx = valsScalar["c"]->get_gradient(q);

  //"n1"fields
  scalarvalueType n1 = valsScalar["n1"]->get_value(q);
  scalargradType n1x = valsScalar["n1"]->get_gradient(q);
  //"n2"fields
  scalarvalueType n2 = valsScalar["n2"]->get_value(q);
  scalargradType n2x = valsScalar["n2"]->get_gradient(q);
  //"n3"fields
  scalarvalueType n3 = valsScalar["n3"]->get_value(q);
  scalargradType n3x = valsScalar["n3"]->get_gradient(q);

  //"u" fields
  vectorgradType ux = valsVector["u"]->get_gradient(q);
  vectorgradType Rux;

  //"u"
  //compute E2=(E-E0)
  vectorgradType E2;
  E2[0][0]=ux[0][0]-(sf1Strain[0][0]*h1V+sf2Strain[0][0]*h2V+sf3Strain[0][0]*h3V);
  E2[1][1]=ux[1][1]-(sf1Strain[1][1]*h1V+sf2Strain[1][1]*h2V+sf3Strain[1][1]*h3V);
  E2[0][1]=constV(0.5)*(ux[0][1]+ux[1][0])-(sf1Strain[0][1]*h1V+sf2Strain[0][1]*h2V+sf3Strain[0][1]*h3V); 
  E2[1][0]=constV(0.5)*(ux[0][1]+ux[1][0])-(sf1Strain[1][0]*h1V+sf2Strain[1][0]*h2V+sf3Strain[1][0]*h3V); 

  //compute stress
  //S=C*(E-E0)
  vectorgradType S;
  computeStress<dim>(CIJ, E2, S);

  //fill residual corresponding to mechanics
  //R=-C*(E-E0)
  for (unsigned int i=0; i<dim; i++){
    for (unsigned int j=0; j<dim; j++){
      Rux[i][j] -= S[i][j]; 
    }
  }

  //"n"
  //compute C*(E-E0)*(Esf*Hn)
  VectorizedArray<double> CEE1=make_vectorized_array(0.0);
  VectorizedArray<double> CEE2=make_vectorized_array(0.0);
  VectorizedArray<double> CEE3=make_vectorized_array(0.0);
  for (unsigned int i=0; i<dim; i++){
    for (unsigned int j=0; j<dim; j++){
      CEE1+=S[i][j]*sf1Strain[i][j];
      CEE2+=S[i][j]*sf2Strain[i][j];
      CEE3+=S[i][j]*sf3Strain[i][j];	
    }
  }
  CEE1*=hn1V;
  CEE2*=hn2V;
  CEE3*=hn3V;
  //compute K*nx
  scalargradType Knx1, Knx2, Knx3;
  for (unsigned int a=0; a<dim; a++) {
    Knx1[a]=0.0;
    Knx2[a]=0.0; 
    Knx3[a]=0.0; 
    for (unsigned int b=0; b<dim; b++){
      Knx1[a]+=Kn1[a][b]*n1x[b];
      Knx2[a]+=Kn2[a][b]*n2x[b];
      Knx3[a]+=Kn3[a][b]*n3x[b];
    }
  }
  
  //check to ensure we are working on the intended field
  //concentration field
  if (this->fields[this->currentFieldIndex].name.compare("c")==0){
    //compute residuals
    valsScalar["c"]->submit_value(rcV,q); 
    valsScalar["c"]->submit_gradient(rcxV,q);
  }
  //n1 fields
  else if (this->fields[this->currentFieldIndex].name.compare("n1")==0){
    //compute residuals
    valsScalar["n1"]->submit_value(rn1V,q);   
    valsScalar["n1"]->submit_gradient(rn1xV,q);
  }
  //n2 fields
  else if (this->fields[this->currentFieldIndex].name.compare("n2")==0){
    //compute residuals
    valsScalar["n2"]->submit_value(rn2V,q);   
    valsScalar["n2"]->submit_gradient(rn2xV,q);
  }
  //n3 fields
  else if (this->fields[this->currentFieldIndex].name.compare("n3")==0){
    //compute residuals
    valsScalar["n3"]->submit_value(rn3V,q);   
    valsScalar["n3"]->submit_gradient(rn3xV,q);
  }
  //displacement field
  else if (this->fields[this->currentFieldIndex].name.compare("u")==0){
    valsVector["u"]->submit_gradient(Rux,q);
  }
}

template <int dim>
void  CoupledCHACMechanicsProblem<dim>::getLHS(typeVector& vals, unsigned int q) const{
  //check to ensure we are working on the intended implicit field
  if (this->fields[this->currentFieldIndex].name.compare("u")==0){
    vectorgradType ux = vals.get_gradient(q);
    vectorgradType Rux;
    
    //compute stress
    vectorgradType S;
    computeStress<dim>(CIJ, ux, S);
    
    //compute residual
    for (unsigned int i=0; i<dim; i++){
      for (unsigned int j=0; j<dim; j++){
	Rux[i][j] = S[i][j]; 
      }
    }
    
    //submit residual value for quadrature integration and assemble
    vals.submit_gradient(Rux,q);
  }
}

//compute integrated free energy value over the domain
template <int dim>
void CoupledCHACMechanicsProblem<dim>::computeFreeEnergyValue(std::vector<double>& freeEnergyValues){
  double value=0.0;
  QGauss<dim>  quadrature_formula(finiteElementDegree+1);
  FE_Q<dim> FE (QGaussLobatto<1>(finiteElementDegree+1));
  FEValues<dim> fe_values (FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);
  const unsigned int   dofs_per_cell = FE.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  std::vector<double> cVal(n_q_points), n1Val(n_q_points), n2Val(n_q_points), n3Val(n_q_points);
  std::vector<Tensor<1,dim,double> > cxVal(n_q_points), n1xVal(n_q_points), n2xVal(n_q_points), n3xVal(n_q_points);

  typename DoFHandler<dim>::active_cell_iterator cell= this->dofHandlersSet[0]->begin_active(), endc = this->dofHandlersSet[0]->end();

  for (; cell!=endc; ++cell) {
	  if (cell->is_locally_owned()){
    	fe_values.reinit (cell);

    	unsigned int fieldIndex;
    	fieldIndex=this->getFieldIndex("c");
    	fe_values.get_function_values(*this->solutionSet[fieldIndex], cVal);
    	fe_values.get_function_gradients(*this->solutionSet[fieldIndex], cxVal);

    	fieldIndex=this->getFieldIndex("n1");
    	fe_values.get_function_values(*this->solutionSet[fieldIndex], n1Val);
    	fe_values.get_function_gradients(*this->solutionSet[fieldIndex], n1xVal);

    	fieldIndex=this->getFieldIndex("n2");
    	fe_values.get_function_values(*this->solutionSet[fieldIndex], n2Val);
    	fe_values.get_function_gradients(*this->solutionSet[fieldIndex], n2xVal);

    	fieldIndex=this->getFieldIndex("n3");
    	fe_values.get_function_values(*this->solutionSet[fieldIndex], n3Val);
    	fe_values.get_function_gradients(*this->solutionSet[fieldIndex], n3xVal);

    	for (unsigned int q=0; q<n_q_points; ++q){
    		double c=cVal[q];
    		double n1 = n1Val[q];
    		double n2 = n2Val[q];
    		double n3 = n3Val[q];


    		double fgrad = 0;
    		for (int i=0; i<dim; i++){
    			for (int j=0; j<dim; j++){
    				fgrad += Kn1[i][j]*n1xVal[q][i]*n1xVal[q][j];
    			}
    		}
    		for (int i=0; i<dim; i++){
    			for (int j=0; j<dim; j++){
    				fgrad += Kn2[i][j]*n2xVal[q][i]*n2xVal[q][j];
    			}
    		}
    		for (int i=0; i<dim; i++){
    			for (int j=0; j<dim; j++){
    				fgrad += Kn3[i][j]*n3xVal[q][i]*n3xVal[q][j];
    			}
    		}
    		fgrad = 0.5*fgrad;
    		double fhomo = (1.0-(h1V+h2V+h3V))*faV + (h1V+h2V+h3V)*fbV;
    		value+=(fhomo+fgrad)*fe_values.JxW(q);
    	}
	  }
  }

  value=Utilities::MPI::sum(value, MPI_COMM_WORLD);
  freeEnergyValues.push_back(value);
}

#endif
