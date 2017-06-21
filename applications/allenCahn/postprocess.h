
template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The order parameter and its derivatives (names here should match those in the macros above)
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > nx = modelVariablesList[0].scalarGrad;

dealii::Tensor<1, dim, dealii::VectorizedArray<double> > pp_field;
pp_field[0] = nx[0];
pp_field[1] = nx[1];


// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
modelResidualsList[0].scalarValueResidual = std::sqrt(pp_field[0]*pp_field[0]+pp_field[1]*pp_field[1]);


}
