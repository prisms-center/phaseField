//implementation to compute stress given the elasticity matrix (CIJ) and strain tensor

#ifndef COMPUTESTRESS_MECHANICS_H
#define COMPUTESTRESS_MECHANICS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "anisotropy.h"

template <int dim>
void computeStress(const dealii::Table<2, double>& CIJ, const dealii::VectorizedArray<double> strain[][dim], dealii::VectorizedArray<double> R[][dim]){
#if problemDIM==3
  dealii::VectorizedArray<double> S[6], E[6];
  E[0]=strain[0][0]; E[1]=strain[1][1]; E[2]=strain[2][2];
  //In Voigt notation: Engineering shear strain=2*strain
  E[3]=strain[1][2]+strain[2][1]; 
  E[4]=strain[0][2]+strain[2][0];
  E[5]=strain[0][1]+strain[1][0];
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
  E[0]=strain[0][0]; E[1]=strain[1][1]; 
  //In Voigt notation: Engineering shear strain=2*strain
  E[2]=strain[0][1]+strain[1][0];
  for (unsigned int i=0; i<3; i++){
    S[i]=0.0;
    for (unsigned int j=0; j<3; j++){
      S[i]+=CIJ(i,j)*E[j];
    }
  }
  R[0][0]=S[0]; R[1][1]=S[1]; 
  R[0][1]=S[2]; R[1][0]=S[2]; 
#elif problemDIM==1
  dealii::Table<1, dealii::VectorizedArray<double> > S[1], E[1];
  E[0]=strain[0][0]; 
  S[0]=CIJ(0,0)*E[0];
  R[0][0]=S[0]; 
#endif			       
}

#endif

