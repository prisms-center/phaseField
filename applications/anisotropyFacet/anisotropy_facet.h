// Number of orientation vectors used to generate anisotropy
#define n_orients 3

template <int dim, int degree>
void
customPDE<dim, degree>::anisotropy(
  const dealii::Tensor<1, dim, dealii::VectorizedArray<double>> &normal,
  dealii::VectorizedArray<double>                               &gamma,
  dealii::Tensor<1, dim, dealii::VectorizedArray<double>>       &dgammadnormal) const
{
  // Orientations
  // Defining orientations in a static array greatly improves performance, but
  // requires specification of n_orients by a macro (as is done here) or by hand

  double orient[n_orients][dim] = {
    {-0.866025403784, -0.5},
    {0.866025403784,  -0.5},
    {0.0,             1.0 }
  };
  // Orientation parameters
  double w[n_orients]     = {50.0, 50.0, 50.0};
  double alpha[n_orients] = {0.6, 0.6, 0.6};

  gamma = constV(1.0);
  for (unsigned int i = 0; i < n_orients; ++i)
    {
      // mn is the dot product of the normal and the orientation vector
      scalarvalueType mn = constV(0.0);
      for (unsigned int j = 0; j < dim; ++j)
        {
          mn += orient[i][j] * normal[j];
        }
      // Application of the heaviside function
      // Vectorized array mn must be unrolled to evaluate conditional
      for (unsigned int j = 0; j < mn.size(); ++j)
        {
          if (mn[j] < 0.0)
            mn[j] = 0.0;
        }
      // Subtracting terms corresponding to the ith orientation from gamma and
      // the components of dgamma/dn
      gamma -= alpha[i] * std::pow(mn, w[i]);
      for (unsigned int j = 0; j < dim; ++j)
        {
          dgammadnormal[j] -= alpha[i] * w[i] * orient[i][j] * std::pow(mn, w[i] - 1.0);
        }
    }
}
