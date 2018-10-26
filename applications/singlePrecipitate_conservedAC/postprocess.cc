// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);

	set_dependencies_value_term_RHS(0, "c, grad(c), n1, grad(n1), grad(u)");
    set_dependencies_gradient_term_RHS(0, "");

    set_output_integral         	(0,true);

	// Variable 1
	set_variable_name				(1,"mu_c");
	set_variable_type				(1,SCALAR);

    set_dependencies_value_term_RHS(1, "c, n1, grad(n1), grad(u)");
    set_dependencies_gradient_term_RHS(1, "");

    set_output_integral         	(1,true);

}

// =================================================================================

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

		scalarvalueType f_chem = (constV(1.0)-(h1V))*faV + (h1V)*fbV + constV(W)*fbarrierV;

		scalarvalueType f_grad = constV(0.0);

		for (int i=0; i<dim; i++){
			for (int j=0; j<dim; j++){
				f_grad += constV(0.5*Kn1[i][j])*n1x[i]*n1x[j];
			}
		}

		// Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
		scalarvalueType cbnV, cbcV, cbcnV,cacV;

		cbcV = faccV/( (constV(1.0)-h1V)*fbccV + h1V*faccV );
		cacV = fbccV/( (constV(1.0)-h1V)*fbccV + h1V*faccV );
		cbnV = hn1V * (c_alpha - c_beta) * cbcV;
		cbcnV = (faccV * (fbccV-faccV) * hn1V)/( ((1.0-h1V)*fbccV + h1V*faccV)*((1.0-h1V)*fbccV + h1V*faccV) );  // Note: this is only true if faV and fbV are quadratic

		// Calculate the stress-free transformation strain and its derivatives at the quadrature point
		dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc, sfts1n, sfts1cn;

		for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
			sfts1[i][j] = constV(sfts_linear1[i][j])*c_beta + constV(sfts_const1[i][j]);
			sfts1c[i][j] = constV(sfts_linear1[i][j]) * cbcV;
			sfts1cc[i][j] = constV(0.0);
			sfts1n[i][j] = constV(sfts_linear1[i][j]) * cbnV;
			sfts1cn[i][j] = constV(sfts_linear1[i][j]) * cbcnV;
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
					CIJ_combined[i][j] = CIJ_Mg[i][j]*(constV(1.0)-h1V) + CIJ_Beta[i][j]*h1V;
				}
			}
			computeStress<dim>(CIJ_combined, E2, S);
		}
		else{
			computeStress<dim>(CIJ_Mg, E2, S);
		}

		scalarvalueType f_el = constV(0.0);

		for (unsigned int i=0; i<dim; i++){
			for (unsigned int j=0; j<dim; j++){
				f_el += constV(0.5) * S[i][j]*E2[i][j];
			}
		}

		total_energy_density = f_chem + f_grad + f_el;

		// Calculate the chemical potential for the concentration
		scalarvalueType mu_c = constV(0.0);
		mu_c += facV*cacV * (constV(1.0)-h1V) + fbcV*cbcV * h1V;
		for (unsigned int i=0; i<dim; i++){
			for (unsigned int j=0; j<dim; j++){
				mu_c -= S[i][j]*( sfts1c[i][j]*h1V);
			}
		}

		// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
		pp_variable_list.set_scalar_value_term_RHS(0, total_energy_density);
		pp_variable_list.set_scalar_value_term_RHS(1, mu_c);

	}
