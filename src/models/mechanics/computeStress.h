// implementation to compute stress given the elasticity matrix (CIJ) and strain
// tensor

#ifndef COMPUTESTRESS_MECHANICS_H
#define COMPUTESTRESS_MECHANICS_H

#include "../../../include/matrixFreePDE.h"

// this source file is temporarily treated as a header file (hence
// #ifndef's) till library packaging scheme is finalized

// #include "anisotropy.h"

// Compute stress with displacement-gradient as input

// Mathematical formulation (3D): S -> stress vector, R -> stress tensor, C ->
// stiffness tensor, E -> strain vector
//
// Generalized Hooke's law:
// [ S(0) ]   [ C11 C12 C13 C14 C15 C16 ] [ E(0) ]
// [ S(1) ]   [ C21 C22 C13 C14 C15 C16 ] [ E(1) ]
// [ S(2) ]   [ C31 C32 C33 C34 C35 C36 ] [ E(2) ]
// [ S(3) ] = [ C41 C42 C43 C44 C45 C46 ] [ E(3) ]
// [ S(4) ]   [ C51 C52 C53 C54 C55 C56 ] [ E(4) ]
// [ S(5) ]   [ C61 C62 C63 C64 C65 C66 ] [ E(5) ]
//
// Strain vector definition: (using Ricci calculus notation, where commas in
// subscript refer to partial derivatives) E(0) = epsilon1 = epsilon11 = 1/2 (
// u1,1 + u1,1 ) E(1) = epsilon2 = epsilon22 = 1/2 ( u2,2 + u2,2 ) E(2) =
// epsilon3 = epsilon33 = 1/2 ( u3,3 + u3,3 ) E(3) = epsilon4 = 2*epsilon23 = (
// u2,3 + u3,2 ) E(4) = epsilon5 = 2*epsilon13 = ( u1,3 + u3,1 ) E(5) = epsilon6
// = 2*epsilon12 = ( u1,2 + u2,1 )
//
// Stress vector definition:
// S(0) = R[0][0] = sigma11
// S(1) = R[1][1] = sigma22
// S(2) = R[2][2] = sigma33
// S(3) = R[1][2] = R[2][1] = sigma23 = sigma32
// S(4) = R[0][2] = R[2][0] = sigma13 = sigma31
// S(5) = R[0][1] = R[1][0] = sigma12 = sigma21

// Currently there are four overloaded versions of this function. They treat the
// cases where CIJ can be a table, a vectorized array, or a tensor and where
// strain and R can be vectorized arrays or tensors. Going forward, there may be
// a better way to reorganize this with templates.

// Overloaded function where CIJ is a table, and the stress and strain are
// vectorized arrays
template <int dim>
void
computeStress(const dealii::Table<2, double>       &CIJ,
              const dealii::VectorizedArray<double> strain[][dim],
              dealii::VectorizedArray<double>       R[][dim])
{
  if (dim == 3)
    {
      dealii::VectorizedArray<double> S[6], E[6];
      E[0] = strain[0][0];
      E[1] = strain[1][1];
      E[2] = strain[2][2];
      // In Voigt notation: Engineering shear strain=2*strain
      E[3] = strain[1][2] + strain[2][1];
      E[4] = strain[0][2] + strain[2][0];
      E[5] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < 6; i++)
        {
          S[i] = 0.0;
          for (unsigned int j = 0; j < 6; j++)
            {
              S[i] += CIJ(i, j) * E[j];
            }
        }
      R[0][0] = S[0];
      R[1][1] = S[1];
      R[2][2] = S[2];
      R[1][2] = S[3];
      R[0][2] = S[4];
      R[0][1] = S[5];
      R[2][1] = S[3];
      R[2][0] = S[4];
      R[1][0] = S[5];

      // Optimized algorithm that skips the zero entries of CIJ for an
      // orthotropic material and is a few percent faster
      //	  dealii::VectorizedArray<double> S[6], E[6];
      //	  E[0]=strain[0][0]; E[1]=strain[1][1]; E[2]=strain[2][2];
      //	  //In Voigt notation: Engineering shear strain=2*strain
      //	  E[3]=strain[1][2]+strain[2][1];
      //	  E[4]=strain[0][2]+strain[2][0];
      //	  E[5]=strain[0][1]+strain[1][0];
      //	  for (unsigned int i=0; i<3; i++){
      //	    S[i]=0.0;
      //	    for (unsigned int j=0; j<3; j++){
      //	      S[i]+=CIJ(i,j)*E[j];
      //	    }
      //	  }
      //
      //	  for (unsigned int i=3; i<6; i++){
      //		  S[i]=CIJ(i,i)*E[i];
      //	  }
      //	  R[0][0]=S[0]; R[1][1]=S[1]; R[2][2]=S[2];
      //	  R[1][2]=S[3]; R[0][2]=S[4]; R[0][1]=S[5];
      //	  R[2][1]=S[3]; R[2][0]=S[4]; R[1][0]=S[5];
    }
  else if (dim == 2)
    {
      dealii::VectorizedArray<double> S[3], E[3];
      E[0] = strain[0][0];
      E[1] = strain[1][1];
      // In Voigt notation: Engineering shear strain=2*strain
      E[2] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < 3; i++)
        {
          S[i] = 0.0;
          for (unsigned int j = 0; j < 3; j++)
            {
              S[i] += CIJ(i, j) * E[j];
            }
        }
      R[0][0] = S[0];
      R[1][1] = S[1];
      R[0][1] = S[2];
      R[1][0] = S[2];
    }
  else
    {
      dealii::VectorizedArray<double> S[1], E[1];
      E[0]    = strain[0][0];
      S[0]    = CIJ(0, 0) * E[0];
      R[0][0] = S[0];
    }
}

// Overloaded function where CIJ, the strain, and the stress are all vectorized
// arrays
template <int dim>
void
computeStress(
  const dealii::VectorizedArray<double> CIJ[2 * dim - 1 + dim / 3][2 * dim - 1 + dim / 3],
  const dealii::VectorizedArray<double> strain[][dim],
  dealii::VectorizedArray<double>       R[][dim])
{
  if (dim == 3)
    {
      dealii::VectorizedArray<double> S[6], E[6];
      E[0] = strain[0][0];
      E[1] = strain[1][1];
      E[2] = strain[2][2];
      // In Voigt notation: Engineering shear strain=2*strain
      E[3] = strain[1][2] + strain[2][1];
      E[4] = strain[0][2] + strain[2][0];
      E[5] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < 6; i++)
        {
          S[i] = 0.0;
          for (unsigned int j = 0; j < 6; j++)
            {
              S[i] += CIJ[i][j] * E[j];
            }
        }
      R[0][0] = S[0];
      R[1][1] = S[1];
      R[2][2] = S[2];
      R[1][2] = S[3];
      R[0][2] = S[4];
      R[0][1] = S[5];
      R[2][1] = S[3];
      R[2][0] = S[4];
      R[1][0] = S[5];
    }
  else if (dim == 2)
    {
      dealii::VectorizedArray<double> S[3], E[3];
      E[0] = strain[0][0];
      E[1] = strain[1][1];
      // In Voigt notation: Engineering shear strain=2*strain
      E[2] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < 3; i++)
        {
          S[i] = 0.0;
          for (unsigned int j = 0; j < 3; j++)
            {
              S[i] += CIJ[i][j] * E[j];
            }
        }
      R[0][0] = S[0];
      R[1][1] = S[1];
      R[0][1] = S[2];
      R[1][0] = S[2];
    }
  else
    {
      dealii::VectorizedArray<double> S[1], E[1];
      E[0]    = strain[0][0];
      S[0]    = CIJ[0][0] * E[0];
      R[0][0] = S[0];
    }
}

// Overloaded function where CIJ is stored as a tensor and the strain and stress
// are vectorized arrays
template <int dim>
void
computeStress(
  const dealii::Tensor<2, 2 * dim - 1 + dim / 3, dealii::VectorizedArray<double>> &CIJ,
  const dealii::VectorizedArray<double> strain[][dim],
  dealii::VectorizedArray<double>       R[][dim])
{
  if (dim == 3)
    {
      dealii::VectorizedArray<double> S[6], E[6];
      E[0] = strain[0][0];
      E[1] = strain[1][1];
      E[2] = strain[2][2];
      // In Voigt notation: Engineering shear strain=2*strain
      E[3] = strain[1][2] + strain[2][1];
      E[4] = strain[0][2] + strain[2][0];
      E[5] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < 6; i++)
        {
          S[i] = 0.0;
          for (unsigned int j = 0; j < 6; j++)
            {
              S[i] += CIJ[i][j] * E[j];
            }
        }
      R[0][0] = S[0];
      R[1][1] = S[1];
      R[2][2] = S[2];
      R[1][2] = S[3];
      R[0][2] = S[4];
      R[0][1] = S[5];
      R[2][1] = S[3];
      R[2][0] = S[4];
      R[1][0] = S[5];
    }
  else if (dim == 2)
    {
      dealii::VectorizedArray<double> S[3], E[3];
      E[0] = strain[0][0];
      E[1] = strain[1][1];
      // In Voigt notation: Engineering shear strain=2*strain
      E[2] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < 3; i++)
        {
          S[i] = 0.0;
          for (unsigned int j = 0; j < 3; j++)
            {
              S[i] += CIJ[i][j] * E[j];
            }
        }
      R[0][0] = S[0];
      R[1][1] = S[1];
      R[0][1] = S[2];
      R[1][0] = S[2];
    }
  else
    {
      dealii::VectorizedArray<double> S[1], E[1];
      E[0]    = strain[0][0];
      S[0]    = CIJ[0][0] * E[0];
      R[0][0] = S[0];
    }
}

// Overloaded function where CIJ, the strain, and the stress are all stored as
// tensors
template <int dim>
void
computeStress(
  const dealii::Tensor<2, 2 * dim - 1 + dim / 3, dealii::VectorizedArray<double>> &CIJ,
  const dealii::Tensor<2, dim, dealii::VectorizedArray<double>>                    strain,
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>>                         &R)
{
  dealii::Tensor<1, 2 * dim - 1 + dim / 3, dealii::VectorizedArray<double>> S, E;

  if (dim == 3)
    {
      E[0] = strain[0][0];
      E[1] = strain[1][1];
      E[2] = strain[2][2];
      // In Voigt notation: Engineering shear strain=2*strain
      E[3] = strain[1][2] + strain[2][1];
      E[4] = strain[0][2] + strain[2][0];
      E[5] = strain[0][1] + strain[1][0];

      // Multiply CIJ and E (in the language of Deal.II this is a tensor
      // contraction) to get S
      S = CIJ * E;

      R[0][0] = S[0];
      R[1][1] = S[1];
      R[2][2] = S[2];
      R[1][2] = S[3];
      R[0][2] = S[4];
      R[0][1] = S[5];
      R[2][1] = S[3];
      R[2][0] = S[4];
      R[1][0] = S[5];
    }
  else if (dim == 2)
    {
      E[0] = strain[0][0];
      E[1] = strain[1][1];
      // In Voigt notation: Engineering shear strain=2*strain
      E[2] = strain[0][1] + strain[1][0];

      // Multiply CIJ and E (in the language of Deal.II this is a tensor
      // contraction) to get S
      S = CIJ * E;

      R[0][0] = S[0];
      R[1][1] = S[1];
      R[0][1] = S[2];
      R[1][0] = S[2];
    }
  else
    {
      R[0][0] = CIJ[0][0] * strain[0][0];
    }
}

// Overloaded function where CIJ is a table and the strain and the stress are
// stored as tensors
template <int dim>
void
computeStress(const dealii::Table<2, double>                                &CIJ,
              const dealii::Tensor<2, dim, dealii::VectorizedArray<double>> &strain,
              dealii::Tensor<2, dim, dealii::VectorizedArray<double>>       &R)
{
  if (dim == 3)
    {
      dealii::VectorizedArray<double> S[6], E[6];
      E[0] = strain[0][0];
      E[1] = strain[1][1];
      E[2] = strain[2][2];
      // In Voigt notation: Engineering shear strain=2*strain
      E[3] = strain[1][2] + strain[2][1];
      E[4] = strain[0][2] + strain[2][0];
      E[5] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < 6; i++)
        {
          S[i] = 0.0;
          for (unsigned int j = 0; j < 6; j++)
            {
              S[i] += CIJ(i, j) * E[j];
            }
        }
      R[0][0] = S[0];
      R[1][1] = S[1];
      R[2][2] = S[2];
      R[1][2] = S[3];
      R[0][2] = S[4];
      R[0][1] = S[5];
      R[2][1] = S[3];
      R[2][0] = S[4];
      R[1][0] = S[5];
    }
  else if (dim == 2)
    {
      dealii::VectorizedArray<double> S[3], E[3];
      E[0] = strain[0][0];
      E[1] = strain[1][1];
      // In Voigt notation: Engineering shear strain=2*strain
      E[2] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < 3; i++)
        {
          S[i] = 0.0;
          for (unsigned int j = 0; j < 3; j++)
            {
              S[i] += CIJ[i][j] * E[j];
            }
        }
      R[0][0] = S[0];
      R[1][1] = S[1];
      R[0][1] = S[2];
      R[1][0] = S[2];
    }
  else
    {
      dealii::VectorizedArray<double> S[1], E[1];
      E[0]    = strain[0][0];
      S[0]    = CIJ[0][0] * E[0];
      R[0][0] = S[0];
    }
}

#endif
