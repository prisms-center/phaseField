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
	double k2;
    double tau;
	if ((p[0] > gbll) && (p[0] < gbrl)){
		k2 = k2_gb;
        tau=tau_gb;
	}
	else {
		k2 = k2_b;
        tau=tau_b;
	}

	double J=k1*exp(-k2/(std::max(ssf,1.0e-6)))*exp(-tau/(this->currentTime));
	double retProb=1.0-exp(-J*userInputs.dtValue*((double)userInputs.steps_between_nucleation_attempts)*dV);
    return retProb;
}
