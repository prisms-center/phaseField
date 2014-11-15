//Compute stress given elasticity matrix, CIJ, and displacement gradient
#ifndef COMPUTESTRESS_H
#define COMPUTESTRESS_H

#include "elasticityModels.h"

template <int dim>
void computeStress(const dealii::Table<2, double>& CIJ, dealii::Tensor<1, dim, gradType>& ux, dealii::VectorizedArray<double> R[][dim]){
#if problemDIM==3
      dealii::VectorizedArray<double> S[6], E[6];
      E[0]=ux[1,0][1,0]; E[1]=ux[1,1][1,1]; E[2]=ux[1,2][1,2];
      E[3]=ux[1,1][1,2]+ux[1,2][1,1];
      E[4]=ux[1,0][1,2]+ux[1,2][1,0];
      E[5]=ux[1,0][1,1]+ux[1,1][1,0];
      for (unsigned int i=0; i<6; i++){
	S[i]=0.0;
	for (unsigned int j=0; j<6; j++){
	  S[i]+=CIJ(i,j)*E[j];
	}
      }
      R[0][0]=S[0]; R[1][1]=S[1]; R[2][2]=S[2];
      R[1][2]=S[3]; R[0][2]=S[4]; R[0][1]=S[5];
      R[2][1]=S[3]; R[2][0]=S[4]; R[1][0]=S[5];     
#elif problemDIM==2
      dealii::VectorizedArray<double> S[3], E[3];
      E[0]=ux[1,0][1,0]; E[1]=ux[1,1][1,1]; 
      E[2]=ux[1,0][1,1]+ux[1,1][1,0];
      for (unsigned int i=0; i<3; i++){
	S[i]=0.0;
	for (unsigned int j=0; j<3; j++){
	  S[i]+=CIJ(i,j)*E[j];
	}
      }
      R[0][0]=S[0]; R[1][1]=S[1]; 
      R[0][1]=S[2]; R[1][0]=S[2]; 
#elif problemDIM==1
      dealii::VectorizedArray<double> S[1], E[1];
      E[0]=ux[1,0]; 
      S[0]=CIJ(0,0)*E[0];
      R[0][0]=S[0];
#endif			       
}


#endif
