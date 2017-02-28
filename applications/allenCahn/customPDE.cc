template <int dim, int degree>
customPDE<dim,degree>::customPDE(userInputParameters _userInputs): MatrixFreePDE<dim,degree>(_userInputs) {


}


template class customPDE<2,1>;
