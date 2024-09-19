// =================================================================================
// NUCLEATION FUNCTIONS
// =================================================================================

// =================================================================================
// Nucleation probability
// =================================================================================
template <int dim, int degree>
double
customPDE<dim, degree>::getNucleationProbability(
  [[maybe_unused]] variableValueContainer variable_value,
  [[maybe_unused]] double                 dV,
  [[maybe_unused]] dealii::Point<dim>     p,
  [[maybe_unused]] unsigned int           variable_index) const
{
  // Supersaturation factor
  double ssf;
  if (dim == 2)
    ssf = variable_value(0) - calmin;
  if (dim == 3)
    ssf = (variable_value(0) - calmin) * (variable_value(0) - calmin);
  // Calculate the nucleation rate
  double J = k1 * exp(-k2 / (std::max(ssf, 1.0e-6))) * exp(-tau / (this->currentTime));
  double retProb =
    1.0 - exp(-J * userInputs.dtValue *
              ((double) userInputs.steps_between_nucleation_attempts) * dV);
  return retProb;
}
