#define n_orients 3
#define gamma_mac (1.0-(0.6*std::pow((-0.866025403784*normal[0]+-0.5*normal[1]),50.)*hs[0]+0.6*std::pow((0.866025403784*normal[0]+-0.5*normal[1]),50.)*hs[1]+0.6*std::pow((0.0*normal[0]+1.0*normal[1]),50.)*hs[2]))
#define gammanx (-(0.6*50.*-0.866025403784*std::pow((-0.866025403784*normal[0]+-0.5*normal[1]),49.)*hs[0]+0.6*50.*0.866025403784*std::pow((0.866025403784*normal[0]+-0.5*normal[1]),49.)*hs[1]+0.6*50.*0.0*std::pow((0.0*normal[0]+1.0*normal[1]),49.)*hs[2]))
#define gammany (-(0.6*50.*-0.5*std::pow((-0.866025403784*normal[0]+-0.5*normal[1]),49.)*hs[0]+0.6*50.*-0.5*std::pow((0.866025403784*normal[0]+-0.5*normal[1]),49.)*hs[1]+0.6*50.*1.0*std::pow((0.0*normal[0]+1.0*normal[1]),49.)*hs[2]))
#define heaviside_macro {(-0.866025403784*normal[0]+-0.5*normal[1]),(0.866025403784*normal[0]+-0.5*normal[1]),(0.0*normal[0]+1.0*normal[1])}

//#include "../../include/typeDefs.h"
//typedef dealii::VectorizedArray<double> scalarvalueType;
//typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double> > scalargradType;

template <int dim,int degree>
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > customPDE<dim,degree>::anisotropy(dealii::Tensor<1, dim, dealii::VectorizedArray<double> > nx) const {

scalarvalueType normgradn = std::sqrt(nx.norm_square());
scalargradType  normal = nx/(normgradn+constV(1.0e-16));
scalarvalueType hs[n_orients] = heaviside_macro;

for (unsigned int i=0; i<n_orients; ++i){
    for (unsigned int j=0; j<hs[i].n_array_elements; ++j){
        if (hs[i][j] < 0.0) hs[i][j] = 0.0;
    }
}

scalarvalueType gamma = gamma_mac;
scalargradType dgammadnormal;
dgammadnormal[0] = gammanx;
dgammadnormal[1] = gammany;

scalargradType aniso;
for (unsigned int i=0; i<dim; ++i){
  for (unsigned int j=0; j<dim; ++j){
      aniso[i] += -normal[i]*normal[j]*dgammadnormal[j];
      if (i==j) aniso[i] +=dgammadnormal[j];
   }
}
aniso = gamma*(aniso*normgradn+gamma*nx);

return aniso;

}
