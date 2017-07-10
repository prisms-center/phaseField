
template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
	variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
	const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

		scalarvalueType total_energy_density = constV(0.0);

		/// The concentration and its derivatives (names here should match those in the macros above)
		scalarvalueType c = variable_list.get_scalar_value(0);
		scalargradType cx = variable_list.get_scalar_gradient(0);

		// The first order parameter and its derivatives (names here should match those in the macros above)
		scalarvalueType n1 = variable_list.get_scalar_value(1);
		scalargradType n1x = variable_list.get_scalar_gradient(1);

		// The derivative of the displacement vector (names here should match those in the macros above)
		vectorgradType ux = variable_list.get_vector_gradient(2);

		scalarvalueType f_chem = (constV(1.0)-(h1V))*faV + (h1V)*fbV;

		scalarvalueType f_grad = constV(0.0);

		for (int i=0; i<dim; i++){
			for (int j=0; j<dim; j++){
				f_grad += constV(0.5*Kn1[i][j])*n1x[i]*n1x[j];
			}
		}

		// Calculate the stress-free transformation strain and its derivatives at the quadrature point
		dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc;

		for (unsigned int i=0; i<dim; i++){
			for (unsigned int j=0; j<dim; j++){
				// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c + b_p
				sfts1[i][j] = constV(sfts_linear1[i][j])*c + constV(sfts_const1[i][j]);
				sfts1c[i][j] = constV(sfts_linear1[i][j]);
				sfts1cc[i][j] = constV(0.0);
			}
		}

		//compute E2=(E-E0)
		dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

		for (unsigned int i=0; i<dim; i++){
			for (unsigned int j=0; j<dim; j++){
				E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V);

			}
		}

		//compute stress
		//S=C*(E-E0)
		dealii::VectorizedArray<double> CIJ_combined[2*dim-1+dim/3][2*dim-1+dim/3];

		if (n_dependent_stiffness == true){
			for (unsigned int i=0; i<2*dim-1+dim/3; i++){
				for (unsigned int j=0; j<2*dim-1+dim/3; j++){
					CIJ_combined[i][j] = userInputs.CIJ_list[0][i][j]*(constV(1.0)-h1V) + userInputs.CIJ_list[1][i][j]*h1V;
				}
			}
			computeStress<dim>(CIJ_combined, E2, S);
		}
		else{
			computeStress<dim>(userInputs.CIJ_list[0], E2, S);
		}

		scalarvalueType f_el = constV(0.0);

		for (unsigned int i=0; i<dim; i++){
			for (unsigned int j=0; j<dim; j++){
				f_el += constV(0.5) * S[i][j]*E2[i][j];
			}
		}

		total_energy_density = f_chem + f_grad + f_el;

		// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
		pp_variable_list.set_scalar_value_residual_term(0, total_energy_density);

	}
