#include "../../include/dealIIheaders.h"

#define problemDIM 3

////define data types

typedef dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > vectorgradType;

void computeStress(const dealii::Table<2, double>& CIJ, vectorgradType& ux, vectorgradType& R);

//template <int dim>
class unitTest
{
public:
	void test_computeStress();
};



