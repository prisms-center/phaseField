// Class for the list of the elasticity tensors, basically just a templated vector of deal.II tensors

#ifndef INCLUDE_LISTOFCIJ_H_
#define INCLUDE_LISTOFCIJ_H_

template<int dim>
class list_of_CIJ
{
public:
	std::vector<dealii::Tensor<2, 2*dim-1+dim/3, dealii::VectorizedArray<double> > > CIJ_list;
};

#endif
