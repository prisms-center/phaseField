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
  void getRHS(const MatrixFree<dim,double> &data, 
	      std::vector<vectorType*> &dst, 
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const;
    
  //LHS implementation for implicit solve 
  void  getLHS(const MatrixFree<dim,double> &data, 
	       vectorType &dst, 
	       const vectorType &src,
	       const std::pair<unsigned int,unsigned int> &cell_range) const;

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
}

template <int dim>
void  CoupledCHACMechanicsProblem<dim>::getRHS(const MatrixFree<dim,double> &data, 
					       std::vector<vectorType*> &dst, 
					       const std::vector<vectorType*> &src,
					       const std::pair<unsigned int,unsigned int> &cell_range) const{
  
  //initialize fields
  typeScalar cVals(data, 0), n1Vals(data,1), n2Vals(data,2), n3Vals(data,3);
  typeVector uVals(data, 4);

  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //initialize c field
    cVals.reinit(cell); cVals.read_dof_values_plain(*src[0]); cVals.evaluate(true, true, false);

    //initialize n fields
    n1Vals.reinit(cell); n1Vals.read_dof_values_plain(*src[1]); n1Vals.evaluate(true, true, false);
    n2Vals.reinit(cell); n2Vals.read_dof_values_plain(*src[2]); n2Vals.evaluate(true, true, false);
    n3Vals.reinit(cell); n3Vals.read_dof_values_plain(*src[3]); n3Vals.evaluate(true, true, false);

    //initialize u field 
    uVals.reinit(cell); uVals.read_dof_values_plain(*src[4]);

    if (c_dependent_misfit == true){
    	uVals.evaluate(false, true, true);
    }
    else{
    	uVals.evaluate(false, true, false);
    }

    //loop over quadrature points
    for (unsigned int q=0; q<cVals.n_q_points; ++q){
      //c
      scalarvalueType c = cVals.get_value(q);
      scalargradType cx = cVals.get_gradient(q);

      //n1
      scalarvalueType n1 = n1Vals.get_value(q);
      scalargradType n1x = n1Vals.get_gradient(q);

      //n2
      scalarvalueType n2 = n2Vals.get_value(q);
      scalargradType n2x = n2Vals.get_gradient(q);

      //n3 
      scalarvalueType n3 = n3Vals.get_value(q);
      scalargradType n3x = n3Vals.get_gradient(q);
      
      //u
      vectorgradType ux = uVals.get_gradient(q);
      vectorgradType Rux;
      
      //if (c_dependent_misfit == true){
    	  vectorhessType uxx = uVals.get_hessian(q);
      //}

      //compute E2=(E-E0)
      dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

      // The result changes if I initialize this to constV(0.0) vs c. That shouldn't be true.
      //double c_effective = 0.125; //c[0];
      dealii::VectorizedArray<double> c_effective = c;

      double cutoff = 0.01;
      if (c[0] < cutoff){
    	  c_effective[0] = cutoff;
      }

      if (c_dependent_misfit == true){
    	  for (unsigned int i=0; i<dim; i++){
    		  for (unsigned int j=0; j<dim; j++){
    			  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( (sf1Strain_linear[i][j]*c_effective+sf1Strain_const[i][j])*h1V + (sf2Strain_linear[i][j]*c_effective+sf2Strain_const[i][j])*h2V + (sf3Strain_linear[i][j]*c_effective+sf3Strain_const[i][j])*h3V);
    		  }
    	  }
      }
      else{
    	  for (unsigned int i=0; i<dim; i++){
    		  for (unsigned int j=0; j<dim; j++){
    			  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-(sf1Strain[i][j]*h1V+sf2Strain[i][j]*h2V+sf3Strain[i][j]*h3V);
    		  }
    	  }
      }
      
      //compute stress
      //S=C*(E-E0)
      computeStress<dim>(CIJ, E2, S);
      
      //fill residual corresponding to mechanics
      //R=-C*(E-E0)
      for (unsigned int i=0; i<dim; i++){
    	  for (unsigned int j=0; j<dim; j++){
    		  Rux[i][j] -= S[i][j];
    	  }
      }
      
      //compute the stress term in the order parameter chemical potential, CEE = C*(E-E0)*(Esf*Hn)
      VectorizedArray<double> CEE1=make_vectorized_array(0.0);
      VectorizedArray<double> CEE2=make_vectorized_array(0.0);
      VectorizedArray<double> CEE3=make_vectorized_array(0.0);
      
      if (c_dependent_misfit == true){
		  for (unsigned int i=0; i<dim; i++){
			  for (unsigned int j=0; j<dim; j++){
				  CEE1+=S[i][j]*(sf1Strain_linear[i][j]*c_effective+sf1Strain_const[i][j]);
				  CEE2+=S[i][j]*(sf2Strain_linear[i][j]*c_effective+sf2Strain_const[i][j]);
				  CEE3+=S[i][j]*(sf3Strain_linear[i][j]*c_effective+sf3Strain_const[i][j]);
			  }
		  }
      }
      else{
    	  for (unsigned int i=0; i<dim; i++){
          	  for (unsigned int j=0; j<dim; j++){
          		  CEE1+=S[i][j]*sf1Strain[i][j];
          		  CEE2+=S[i][j]*sf2Strain[i][j];
          		  CEE3+=S[i][j]*sf3Strain[i][j];
          	  }
            }
      }
      CEE1*=hn1V;
      CEE2*=hn2V;
      CEE3*=hn3V;
      
      // compute the stress term in the concentration chemical potential, CEEcx = [C*(E-E0)*E0c]x, must be a vector with length dim
      dealii::VectorizedArray<double> CEEcx[dim];

      for (unsigned int k=0; k<dim; k++){
    	  CEEcx[k] = constV(0.0);
      }

      if (c[0] >= cutoff){
      if (c_dependent_misfit == true){
    	  dealii::VectorizedArray<double> E3[dim][dim], S3[dim][dim];

    		  for (unsigned int i=0; i<dim; i++){
    			  for (unsigned int j=0; j<dim; j++){
    				  E3[i][j] =  (sf1Strain_linear[i][j]*h1V + sf2Strain_linear[i][j]*h2V + sf3Strain_linear[i][j]*h3V);
    			  }
    		  }

    		  computeStress<dim>(CIJ, E3, S3);
    		  for (unsigned int k=0; k<dim; k++){
    			  for (unsigned int i=0; i<dim; i++){
    				  for (unsigned int j=0; j<dim; j++){
    					  CEEcx[k]+=S3[i][j] * (constV(0.5)*(uxx[i][j][k]+uxx[j][i][k]) - (sf1Strain_linear[i][j]*h1V + sf2Strain_linear[i][j]*h2V + sf3Strain_linear[i][j]*h3V)*cx[k]);
    				  }
    			  }
    		  }
      }
      }

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
  
      //submit values
      cVals.submit_value(rcV,q); cVals.submit_gradient(rcxV,q);
      n1Vals.submit_value(rn1V,q); n1Vals.submit_gradient(rn1xV,q);
      n2Vals.submit_value(rn2V,q); n2Vals.submit_gradient(rn2xV,q);
      n3Vals.submit_value(rn3V,q); n3Vals.submit_gradient(rn3xV,q);
      uVals.submit_gradient(Rux,q);
    }
    
    //integrate values
    cVals.integrate(true, true);  cVals.distribute_local_to_global(*dst[0]);
    n1Vals.integrate(true, true); n1Vals.distribute_local_to_global(*dst[1]);
    n2Vals.integrate(true, true); n2Vals.distribute_local_to_global(*dst[2]);
    n3Vals.integrate(true, true); n3Vals.distribute_local_to_global(*dst[3]);
    uVals.integrate(false, true); uVals.distribute_local_to_global(*dst[4]);
  }
}

template <int dim>
void  CoupledCHACMechanicsProblem<dim>::getLHS(const MatrixFree<dim,double> &data, 
					       vectorType &dst, 
					       const vectorType &src,
					       const std::pair<unsigned int,unsigned int> &cell_range) const{
  typeVector uVals(data, 4);
  
  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //initialize u field 
    uVals.reinit(cell); uVals.read_dof_values_plain(src); uVals.evaluate(false, true, false);

    //loop over quadrature points
    for (unsigned int q=0; q<uVals.n_q_points; ++q){
      //u
      vectorgradType ux = uVals.get_gradient(q);
      vectorgradType Rux;

      //compute strain tensor
      dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);
	}
      }
    
      //compute stress tensor
      computeStress<dim>(CIJ, E, S);
      
      //compute residual
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  Rux[i][j] = S[i][j]; 
	}
      }
      
      //submit residual value
      uVals.submit_gradient(Rux,q);
    }
    
    //integrate
    uVals.integrate(false, true); uVals.distribute_local_to_global(dst);
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
  std::vector<std::vector<Tensor<1,dim,double> > > uxVal(n_q_points,std::vector<Tensor<1,dim,double> >(dim)); // I would think that this should be a rank 2 tensor, but that doesn't compile
  //std::vector<std::vector<Tensor<1,dim,double> > > uxVal[n_q_points][dim];
  //std::vector<vectorgradType> uxVal(n_q_points);

  // remove later
  double value_homo = 0, value_grad = 0, value_el = 0;

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

    	//fieldIndex=this->getFieldIndex("u");
    	//fe_values.get_function_gradients(*this->solutionSet[fieldIndex], uxVal); // This step doesn't work, uxVal not the right type

    	for (unsigned int q=0; q<n_q_points; ++q){
    		double c=cVal[q];
    		double n1 = n1Val[q];
    		double n2 = n2Val[q];
    		double n3 = n3Val[q];

//    		vectorgradType ux;
//    		for (unsigned int i=0; i<dim; i++){
//    			for (unsigned int j=0; j<dim; j++){
//    				ux[i][j] = uxVal[q][i][j];
//    			}
//    		}

    		// calculate the interfacial energy
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
    		fgrad = 0.5*fgrad + W*fbarrierV; // need to generalize for multiple order parameters

    		// calculate the homogenous chemical energy
    		double fhomo = (1.0-(h1V+h2V+h3V))*faV + (h1V+h2V+h3V)*fbV;

    		// Calculatate the elastic energy
    		double fel = 0;

//    		VectorizedArray<double> test;
//    		test = constV(1.5);
//    		double test_2 = test[0];

/*
    		vectorgradType E2;
    		//dealii::Table<2, double> E2;
    		for (unsigned int i=0; i<dim; i++){
    			for (unsigned int j=0; j<dim; j++){
    				E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-(sf1Strain[i][j]*h1V+sf2Strain[i][j]*h2V+sf3Strain[i][j]*h3V);
    			}
    		}

			#if problemDIM==3
    		dealii::Table<1, dealii::VectorizedArray<double> > E(6);

    		E(0)=E2[0][0]; E(1)=E2[1][1]; E(2)=E2[2][2];
    		E(3)=E2[1][2]+E2[2][1];
    		E(4)=E2[0][2]+E2[2][0];
    		E(5)=E2[0][1]+E2[1][0];

    		for (unsigned int i=0; i<6; i++){
    			for (unsigned int j=0; j<6; j++){
    				fel+=CIJ(i,j)*E(i)[0]*E(j)[0];
    			}
    		}
			#elif problemDIM==2
			  dealii::Table<1, dealii::VectorizedArray<double> > E(3);
			  E(0)=E2[0][0]; E(1)=E2[1][1];
			  E(2)=E2[0][1]+E2[1][0];

			  for (unsigned int i=0; i<3; i++){
				for (unsigned int j=0; j<3; j++){
					fel+=CIJ(i,j)*E(i).operator[](0)*E(j).operator[](0);
				}
			  }
			#elif problemDIM==1
			  dealii::Table<1, dealii::VectorizedArray<double> > E(1);
			  E(0)=E2[0][0];
			  fel=CIJ(0,0)*E(0).operator[](0);
			#endif
			*/


    		// Sum the energies at each integration point
    		value+=(fhomo+fgrad+fel)*fe_values.JxW(q);

    		// Remove this later, sum the components of the energies at each integration point
    		value_homo += fhomo*fe_values.JxW(q);
    		value_grad += fgrad*fe_values.JxW(q);
    		value_el += fel*fe_values.JxW(q);
    	}
	  }
  }

  value=Utilities::MPI::sum(value, MPI_COMM_WORLD);
  //freeEnergyValues.push_back(value);

  // remove later
  value_homo=Utilities::MPI::sum(value_homo, MPI_COMM_WORLD);
  value_grad=Utilities::MPI::sum(value_grad, MPI_COMM_WORLD);
  value_el=Utilities::MPI::sum(value_el, MPI_COMM_WORLD);

  freeEnergyValues.push_back(value_grad);

  std::cout<<"Homogenous Free Energy: "<<value_homo<<std::endl;
  std::cout<<"Interfacial Free Energy: "<<value_grad<<std::endl;
  std::cout<<"Elastic Free Energy: "<<value_el<<std::endl;
}

//structure representing each nucleus
struct nucleus{
  unsigned int index;
  dealii::Point<problemDIM> center;
  double radius;
  double seededTime, seedingTime;
};
//vector of all nucleus seeded in the problem
std::vector<nucleus> nuclei, localNuclei;

//nucleation model implementation
template <int dim>
void CoupledCHACMechanicsProblem<dim>::modifySolutionFields()
{
  //current time
  double t=this->currentTime;
  unsigned int inc=this->currentIncrement;
  double dx=spanX/( (double)subdivisionsX )/std::pow(2.0,refineFactor);
  double rand_val;
  int count = 0;
  //nucleation parameters
  double nRadius = 2.5; //spanX/20.0;
  double minDistBetwenNuclei=4*nRadius;

  unsigned int maxNumberNuclei=5; // doesn't do anything currently

  //get the list of node points in the domain
  std::map<dealii::types::global_dof_index, dealii::Point<dim> > support_points;
  dealii::DoFTools::map_dofs_to_support_points (dealii::MappingQ1<dim>(), *this->dofHandlersSet[0], support_points);
  //fields
  vectorType* n1=this->solutionSet[this->getFieldIndex("n1")];
  vectorType* n2=this->solutionSet[this->getFieldIndex("n2")];
  vectorType* n3=this->solutionSet[this->getFieldIndex("n3")];
  vectorType* c=this->solutionSet[this->getFieldIndex("c")];
  const double k1 = 0.0001; // nucleation probability constant
  const double k2 = 1.0;	// nucleation probability constant
  const double c0 = 0.300;	// baseline concentration?
  double J = 0.0;
  //delete the previous entries in the nuclei vector (old nucleus are still retained in the localNuclei vector)
  nuclei.clear();

  //populate localNuclei vector
  if (inc <= timeIncrements){
    nucleus* temp;
    //add nuclei based on concentration field values
    //loop over all points in the domain
    for (typename std::map<dealii::types::global_dof_index, dealii::Point<dim> >::iterator it=support_points.begin(); it!=support_points.end(); ++it){
      unsigned int dof=it->first;
      //set only local owned values of the parallel vector (eventually turn this into a separate function for each order parameter)
      if (n1->locally_owned_elements().is_element(dof)){
    	  dealii::Point<dim> nodePoint=it->second;
    	  double n1Value=(*n1)(dof);
    	  double n2Value=(*n2)(dof);
    	  double n3Value=(*n3)(dof);
    	  double cValue=(*c)(dof);

    	  rand_val = (double)rand()/(RAND_MAX);

    	  if ((t > 1000000000*timeStep) || (n1Value+n2Value+n3Value > 1.0e-6) || (cValue <= 0.0)) {
    		  J = 0;
    	  }
		  else{
			  J = cValue/avg_Nd * dx*dx/((double)spanX * (double)spanY) * 0.01;
    	  }

    	  if (rand_val <= J){
    		  bool isClose=false;
    		  for (std::vector<nucleus>::iterator thisNuclei=localNuclei.begin(); thisNuclei!=localNuclei.end(); ++thisNuclei){
    			  if (thisNuclei->center.distance(nodePoint)<minDistBetwenNuclei){
    				  isClose=true;
    			  }
    		  }

    		  if (!isClose){
    			  temp = new nucleus;
    			  temp->index=localNuclei.size();
    			  temp->center=nodePoint;
    			  temp->radius=nRadius;
    			  temp->seededTime=t;
    			  temp->seedingTime = 10000.0*timeStep;
    			  localNuclei.push_back(*temp);
    		  }
    	  }
      }
    }


    //filter nuclei by comparing with other processors
    int numProcs=Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    int thisProc=Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    std::vector<int> numNucleiInProcs(numProcs, 0);
    //send nuclei information to processor 0
    int numNuclei=localNuclei.size();
    //send information about number of nuclei to processor 0
    if (thisProc!=0){
      MPI_Send(&numNuclei, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    else{
      numNucleiInProcs[0]=numNuclei;
      for (int proc=1; proc<numProcs; proc++){
	MPI_Recv(&numNucleiInProcs[proc], 1, MPI_INT, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //filter nuclei in processor zero
    //receive nuclei info from all processors
    if (thisProc!=0){
    	if (numNuclei>0){
    		std::vector<double> tempData((dim+3)*numNuclei);
    		unsigned int i=0;
    		for (std::vector<nucleus>::iterator thisNuclei=localNuclei.begin(); thisNuclei!=localNuclei.end(); ++thisNuclei){
    			tempData[i*(dim+3)]=thisNuclei->radius;
    			tempData[i*(dim+3)+1]=thisNuclei->seededTime;
    			tempData[i*(dim+3)+2]=thisNuclei->seedingTime;
    			for (unsigned int j=0; j<dim; j++) tempData[i*(dim+3)+3+j]=thisNuclei->center[j];
    			i++;
    		}
    		MPI_Send(&tempData[0], numNuclei*(dim+3), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    	}
    }
    else{
    	//temporary array to store all the nuclei
    	std::vector<std::vector<double>*> tempNuceli(numProcs);
    	for (int proc=0; proc<numProcs; proc++) {
    		std::vector<double>* temp=new std::vector<double>(numNucleiInProcs[proc]*(dim+3));
    		if (numNucleiInProcs[proc]>0){
    			if (proc==0){
    				unsigned int i=0;
    				for (std::vector<nucleus>::iterator thisNuclei=localNuclei.begin(); thisNuclei!=localNuclei.end(); ++thisNuclei){
    					(*temp)[i*(dim+3)]=thisNuclei->radius;
    					(*temp)[i*(dim+3)+1]=thisNuclei->seededTime;
    					(*temp)[i*(dim+3)+2]=thisNuclei->seedingTime;
    					for (unsigned int j=0; j<dim; j++) (*temp)[i*(dim+3)+3+j]=thisNuclei->center[j];
    					i++;
    				}
    			}
    			else{
    				MPI_Recv(&((*temp)[0]), numNucleiInProcs[proc]*(dim+3), MPI_DOUBLE, proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    			}
    			tempNuceli[proc]=temp;
    		}
    	}

    	//filter the nuclei and add to nuclei vector in processor zero
    	for (int proc1=0; proc1<numProcs; proc1++) {
    		for (int i1=0; i1<numNucleiInProcs[proc1]; i1++){
    			double rad1=(*tempNuceli[proc1])[i1*(dim+3)];
    			double time1=(*tempNuceli[proc1])[i1*(dim+3)+1];
    			double seedingTime1=(*tempNuceli[proc1])[i1*(dim+3)+2];
    			dealii::Point<dim> center1;
    			for (unsigned int j1=0; j1<dim; j1++) {
    				center1[j1]=(*tempNuceli[proc1])[i1*(dim+3)+3+j1];
    			}
    			bool addNuclei=true;
    			//check if this nuceli present in any other processor
    			for (int proc2=0; proc2<numProcs; proc2++) {
    				if (proc1!=proc2){
    					for (int i2=0; i2<numNucleiInProcs[proc2]; i2++){
    						double rad2=(*tempNuceli[proc2])[i2*(dim+3)];
    						double time2=(*tempNuceli[proc2])[i2*(dim+3)+1];
    						dealii::Point<dim> center2;
    						for (unsigned int j2=0; j2<dim; j2++) center2(j2)=(*tempNuceli[proc2])[i2*(dim+3)+3+j2];
    						if ((center1.distance(center2)<=minDistBetwenNuclei) && (time1>=time2)){
    							addNuclei=false;
    							break;
    						}
    					}
    					if (!addNuclei) {break;}
    				}
    			}
    			if (addNuclei){
    				temp = new nucleus;
    				temp->index=nuclei.size();
    				temp->radius=rad1;
    				temp->seededTime=time1;
    				temp->seedingTime=seedingTime1;
    				temp->center=center1;
    				nuclei.push_back(*temp);
    			}
    		}
    	}
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //disperse nuclei to all other processors
    unsigned int numGlobalNuclei;
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0) {numGlobalNuclei=nuclei.size();}
    MPI_Bcast(&numGlobalNuclei, 1, MPI_INT, 0, MPI_COMM_WORLD);
    this->pcout << "total number of nuclei currently seeded : "  << numGlobalNuclei << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    //
    std::vector<double> temp2(numGlobalNuclei*(dim+3));
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0){
      unsigned int i=0;
      for (std::vector<nucleus>::iterator thisNuclei=nuclei.begin(); thisNuclei!=nuclei.end(); ++thisNuclei){
	temp2[i*(dim+3)]=thisNuclei->radius;
	temp2[i*(dim+3)+1]=thisNuclei->seededTime;
	temp2[i*(dim+3)+2]=thisNuclei->seedingTime;
	for (unsigned int j=0; j<dim; j++) temp2[i*(dim+3)+3+j]=thisNuclei->center[j];
	i++;
      }
    }
    MPI_Bcast(&temp2[0], numGlobalNuclei*(dim+3), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //receive all nuclei
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)!=0){
    	for(unsigned int i=0; i<numGlobalNuclei; i++){
    		temp = new nucleus;
    		temp->index=nuclei.size();
    		temp->radius=temp2[i*(dim+3)];
    		temp->seededTime=temp2[i*(dim+3)+1];
    		temp->seedingTime=temp2[i*(dim+3)+2];
    		dealii::Point<dim> tempCenter;
    		for (unsigned int j=0; j<dim; j++) tempCenter[j]=temp2[i*(dim+3)+3+j];
    		temp->center=tempCenter;
    		nuclei.push_back(*temp);
      }
    }
  }

  //seed nuclei
  unsigned int fieldIndex=this->getFieldIndex("n1");
  for (std::vector<nucleus>::iterator thisNuclei=nuclei.begin(); thisNuclei!=nuclei.end(); ++thisNuclei){

	  dealii::Point<dim> center=thisNuclei->center;
	  double radius=thisNuclei->radius;
	  double seededTime=thisNuclei->seededTime;
	  double seedingTime=thisNuclei->seedingTime;
	  this->pcout << "times: " << t << " " << seededTime << " " << seedingTime << std::endl;
	  //loop over all points in the domain
	  for (typename std::map<dealii::types::global_dof_index, dealii::Point<dim> >::iterator it=support_points.begin(); it!=support_points.end(); ++it){
		  unsigned int dof=it->first;
		  //set only local owned values of the parallel vector
		  if (n1->locally_owned_elements().is_element(dof)){
			  dealii::Point<dim> nodePoint=it->second;
			  //check conditions and seed nuclei
			  double r=nodePoint.distance(center);
			  if (r<=(2*radius)){
				  if ((t>seededTime) && (t<(seededTime+seedingTime))){
					  //this->pcout << "times: " << t << " " << seededTime << " " << seedingTime << std::endl;
					  //(*n1)(dof)=0.5*(1.0-std::tanh((r-radius)/(dx)));
					  (*n1)(dof)=0.5*(1.0-std::tanh((r-radius)/(0.4)));
				  }
			  }
		  }
	  }
  }
}

#endif
