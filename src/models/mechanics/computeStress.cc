//implementation to compute stress given elasticity matrix (CIJ) and
//strain/displacement-gradient

#ifndef COMPUTESTRESS_MECHANICS_H
#define COMPUTESTRESS_MECHANICS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "anisotropy.cc"

//Compute stress with displacement-gradient as input
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

//Compute stress with strain as input
template <int dim>
void computeStress(const dealii::Table<2, double>& CIJ, const dealii::VectorizedArray<double> E2[][dim], dealii::VectorizedArray<double> R[][dim]){
#if problemDIM==3
  dealii::VectorizedArray<double> S[6], E[6];
  E[0]=E2[0][0]; E[1]=E2[1][1]; E[2]=E2[2][2];
  E[3]=E2[1][2]; E[4]=E2[0][2]; E[5]=E2[0][1];
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
  E[0]=E2[0][0]; E[1]=E2[1][1]; E[2]=E2[0][1];
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
  E[0]=E2[0,0]; 
  S[0]=CIJ(0,0)*E[0];
  R[0][0]=S[0];
#endif			       
}

#endif
