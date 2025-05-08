// Nucleation functions: nucleation probability

template<int dim, int degree>
double
customPDE<dim, degree>::getNucleationProbability(variableValueContainer variable_value,
    double                 dV,
    Point<dim>             p,
    unsigned int variable_index) const
{
  // TODO (Phil): enable submission of the following several values,
  //               via user input parameters

  // Number of order parameters
  unsigned int number_ops_ini = 23;

  // Number of empty OPs (for nucleation)
  unsigned int number_ops_nuc = 60;

  // Number of dislocation density fields
  unsigned int N_rho = 83;

  // Total number of order parameters
  unsigned int number_ops_total = number_ops_ini + number_ops_nuc;

  double nuclProb = 0.0;

  double rho = variable_value(number_ops_total + N_rho);
  double mult2Op = variable_value(number_ops_total + N_rho + 1);

  double Prob_rho, Prob_gb;

  // The big number in the numerator coming frmo Rayleigh nucleation model and plugging in physical parameters for Mg, the calculation uses SI units
  Prob_rho = std::exp(-(M_PI / 4) * std::pow(33790 / (A * sqrt(rho * 1.0e6)), 2));

  Prob_gb = mult2Op;

  if (Prob_gb > 0.15)
    {
      nuclProb = Prob_rho;
    }

  return nuclProb;
}
