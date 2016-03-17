#include "../../include/dealIIheaders.h"
#include <iostream>

#define finiteElementDegree 1
#define vectorgradType dealii::Tensor<2, dim, dealii::VectorizedArray<double> >
#define typeScalar dealii::FEEvaluation<dim,finiteElementDegree,finiteElementDegree+1,1,double>
#define typeVector dealii::FEEvaluation<dim,finiteElementDegree,finiteElementDegree+1,dim,double>

//define data type

template <int dim>
void computeStress(const dealii::Table<2, double>& CIJ, const dealii::VectorizedArray<double> ux[][dim], const dealii::VectorizedArray<double> R[][dim]);

#include "../../src/models/mechanics/computeStress.h"

template <int dim>
class unitTest
{
public:
	bool test_computeStress();
	bool test_getRHS();
	bool test_computeRHS();
};

#include "test_computeStress.h"
#include "test_getRHS.h"
//#include "test_computeRHS.h"
