// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);

	set_need_value_residual_term	(0,true);
	set_need_gradient_residual_term	(0,false);

	set_output_integral         	(0,true);



	// Variable 2
	set_variable_name				(1,"von_mises_stress");
	set_variable_type				(1,SCALAR);

	set_need_value_residual_term	(1,true);
	set_need_gradient_residual_term	(1,false);
	set_output_integral         	(1,false);




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
		scalarvalueType n1 = variable_list.get_scalar_value(2);
		scalargradType n1x = variable_list.get_scalar_gradient(2);

		// The derivative of the displacement vector (names here should match those in the macros above)
		vectorgradType ux = variable_list.get_vector_gradient(3);

		scalarvalueType faV = (A0+A1*c+A2*c*c+A3*c*c*c+A4*c*c*c*c);
		scalarvalueType facV = (A1 + 2.0*A2*c + 3.0*A3*c*c + 4.0*A4*c*c*c);
		scalarvalueType faccV = (2.0*A2 + 6.0*c + 12.0*A4*c*c);
		scalarvalueType fbV = (B2*c*c + B1*c + B0);
		scalarvalueType fbcV = (2.0*B2*c + B1);
		scalarvalueType fbccV = constV(2.0*B2);
		scalarvalueType h1V = 3.0*n1*n1-2.0*n1*n1*n1; //(10.0*n1*n1*n1-15.0*n1*n1*n1*n1+6.0*n1*n1*n1*n1*n1);
		scalarvalueType hn1V = 6.0*n1-6.0*n1*n1; //(30.0*n1*n1-60.0*n1*n1*n1+30.0*n1*n1*n1*n1);

		scalarvalueType fbarrierV = n1*n1*(1.0 - 2.0*n1 + n1*n1);


		scalarvalueType f_chem = (constV(1.0)-(h1V))*faV + (h1V)*fbV + constV(W)*fbarrierV;

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


		// The Von Mises Stress
		dealii::VectorizedArray<double> vm_stress;
		if (dim == 3){
			vm_stress = (S[0][0]-S[1][1])*(S[0][0]-S[1][1]) + (S[1][1]-S[2][2])*(S[1][1]-S[2][2]) + (S[2][2]-S[0][0])*(S[2][2]-S[0][0]);
			vm_stress += constV(6.0)*(S[0][1]*S[0][1] + S[1][2]*S[1][2] + S[2][0]*S[2][0]);
			vm_stress *= constV(0.5);
			vm_stress = std::sqrt(vm_stress);
		}
		else {
			vm_stress = S[0][0]*S[0][0] - S[0][0]*S[1][1] + S[1][1]*S[1][1] + constV(3.0)*S[0][1]*S[0][1];
			vm_stress = std::sqrt(vm_stress);
		}



		// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
		pp_variable_list.set_scalar_value_residual_term(0, total_energy_density);
		pp_variable_list.set_scalar_value_residual_term(1, vm_stress);

	}
