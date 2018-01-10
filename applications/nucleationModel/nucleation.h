// =================================================================================
// NUCLEATION FUNCTIONS
// =================================================================================

// =================================================================================
// Nucleation probability
// =================================================================================
template <int dim, int degree>
double customPDE<dim,degree>::getNucleationProbability(variableValueContainer variable_value, double dV, dealii::Point<dim> p) const
{
	//Supersaturation factor
    double ssf;
    if (dim ==2) ssf=variable_value(0)-calmin;
    if (dim ==3) ssf=(variable_value(0)-calmin)*(variable_value(0)-calmin);
	// Calculate the nucleation rate
	double J=k1*exp(-k2/(std::max(ssf,1.0e-6)))*exp(-tau/(this->currentTime));
	double retProb=1.0-exp(-J*userInputs.dtValue*((double)userInputs.nucleation_parameters_list[0].steps_between_attempts)*dV);
    return retProb;
}
