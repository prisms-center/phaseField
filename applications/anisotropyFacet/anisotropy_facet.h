// Number of orientation vectors used to generate anisotropy
#define n_orients 3

template <int dim,int degree>
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > customPDE<dim,degree>::anisotropy(dealii::Tensor<1, dim, dealii::VectorizedArray<double> > nx) const {

// Orientations
// Defining orientations in a static array greatly improves performance, but requires
// specification of n_orients by a macro (as is done here) or by hand
double orient[n_orients][dim] = {{-0.866025403784,-0.5},{0.866025403784,-0.5},{0.0,1.0}};
// Orientation parameters
double w[n_orients] = {50.0,50.0,50.0};
double alpha[n_orients] = {0.6,0.6,0.6};
// Calculation of normal vector
scalarvalueType normgradn = std::sqrt(nx.norm_square());
scalargradType  normal = nx/(normgradn+constV(1.0e-16));

scalarvalueType gamma = constV(1.0);
scalargradType dgammadnormal; // this is automatically initialized to (0.0,0.0,0.0)

for (unsigned int i=0; i<n_orients; ++i){
// mn is the dot product of the normal and the orientation vector
    scalarvalueType mn = constV(0.0);
    for (unsigned int j=0; j<dim; ++j){
        mn += orient[i][j]*normal[j];
    }
// Application of the heaviside function
// Vectorized array mn must be unrolled to evaluate conditional
    for (unsigned int j=0; j<mn.n_array_elements; ++j){
        if (mn[j] < 0.0) mn[j] = 0.0;
    }
// Subtracting terms corresponding to the ith orientation from gamma and the
// components of dgamma/dn
    gamma -= alpha[i]*std::pow(mn,w[i]);
    for (unsigned int j=0; j<dim; ++j){
        dgammadnormal[j] -= alpha[i]*w[i]*orient[i][j]*std::pow(mn,w[i]-1.0);
    }
}

// Calculation of anisotropic gradient, exactly as in CHAC_anisotropy
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
