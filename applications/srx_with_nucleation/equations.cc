// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for
// each function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void
customAttributeLoader::loadVariableAttributes()
{
  // TODO: enable submission/determination of the following several values
  //       via user input parameters

  // Number of order parameters
  unsigned int number_ops_ini = 23;

  // Number of empty OPs (for nucleation)
  unsigned int number_ops_nuc = 60;

  // Number of dislocation density fields
  unsigned int N_rho = 83;

  // Total number of order parameters
  unsigned int number_ops_total = number_ops_ini + number_ops_nuc;

  // Fields for initial grains
  for (unsigned int var_index = 0; var_index < number_ops_ini; var_index++)
    {
      std::string var_name = "n";
      var_name.append(std::to_string(var_index));
      std::string grad_var_name = "grad(" + var_name + ")";

      set_variable_name(var_index, var_name);
      set_variable_type(var_index, SCALAR);
      set_variable_equation_type(var_index, EXPLICIT_TIME_DEPENDENT);

      set_dependencies_value_term_RHS(var_index, var_name);
      set_dependencies_gradient_term_RHS(var_index, grad_var_name);
    }

  // Fields for nucleating grains
  for (unsigned int var_index = number_ops_ini; 
       var_index < number_ops_total; var_index++)
    {
      std::string var_name = "n";
      var_name.append(std::to_string(var_index));
      std::string grad_var_name = "grad(" + var_name + ")";

      set_variable_name(var_index, var_name);
      set_variable_type(var_index, SCALAR);
      set_variable_equation_type(var_index, EXPLICIT_TIME_DEPENDENT);

      set_dependencies_value_term_RHS(var_index, var_name);
      set_dependencies_gradient_term_RHS(var_index, grad_var_name);

      set_allowed_to_nucleate(var_index, true);
      set_need_value_nucleation(var_index, false);
    }

  // Fields for dislocation density
  for (unsigned int rho_index = 0; rho_index < N_rho; rho_index++)
    {
      unsigned int var_index = rho_index + number_ops_total;

      std::string var_name = "rho";
      var_name.append(std::to_string(rho_index));

      set_variable_name(var_index, var_name);
      set_variable_type(var_index, SCALAR);
      set_variable_equation_type(var_index, EXPLICIT_TIME_DEPENDENT);

      set_dependencies_value_term_RHS(var_index, var_name);
    }

  // Other variables to track
  unsigned int last_rho_index = number_ops_total + N_rho - 1;

  // Overall dislocation density, interpolated
  set_variable_name(last_rho_index + 1, "rho");
  set_variable_type(last_rho_index + 1, SCALAR);
  set_variable_equation_type(last_rho_index + 1, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(last_rho_index + 1, "rho");
  set_allowed_to_nucleate(last_rho_index + 1, false);
  set_need_value_nucleation(last_rho_index + 1, true);

  // Sum of multiplication of 2 order parameters (double product)
  set_variable_name(last_rho_index + 2, "mult2Op");
  set_variable_type(last_rho_index + 2, SCALAR);
  set_variable_equation_type(last_rho_index + 2, AUXILIARY);

  set_dependencies_value_term_RHS(
    last_rho_index + 2,
    "n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, "
    "n18, n19, n20, n21, n22, n23, multi2Op");

  set_allowed_to_nucleate(last_rho_index + 2, false);
  set_need_value_nucleation(last_rho_index + 2, true);
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time
// dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a
// list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one
// proportional to the test function and one proportional to the gradient of the
// test function. The index for each variable in this list corresponds to the
// index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // --- Getting the values and derivatives of the model variables ---

  // Number of order parameters
  unsigned int number_ops_ini = 23;

  // Number of empty OPs (for nucleation)
  unsigned int number_ops_nuc = 60;

  // Number of dislocation density fields
  unsigned int N_rho = 83;

  // Total number of order parameters
  unsigned int number_ops_total = number_ops_ini + number_ops_nuc;

  VectorizedArray<double> fnV = constV(0.0);
  VectorizedArray<double> fsV = constV(0.0);
  VectorizedArray<double> sum_nsq = constV(0.0);
  
  scalarvalueType rhocalc = constV(0.0);
  scalarvalueType ni, nj, rhoi;
  scalargradType nix;

  std::vector<scalarvalueType> value_terms;
  value_terms.resize(number_ops_total);

  std::vector<scalargradType> gradient_terms;
  gradient_terms.resize(number_ops_total);

  // Nucleation expressions
  std::vector<scalarvalueType> source_term;
  source_term.resize(number_ops_total);
  std::vector<scalarvalueType> gamma;
  gamma.resize(number_ops_total);

  // Initializing source and gamma terms
  for (unsigned int i = 0; i < number_ops_total; i++)
    {
      source_term[i] = constV(0.0);
      gamma[i] = constV(1.0);
    }
  seedNucleus(q_point_loc, source_term, gamma);

  // Calculating rho
  // Sum of squares of order parameters
  for (unsigned int i = 0; i < number_ops_total; i++)
    {
      ni = variable_list.get_scalar_value(i);
      // Capping value of ni to make sure it is not negative
      ni = std::max(ni, constV(0.0));
      sum_nsq += ni * ni;
    }

  for (unsigned int i = 0; i < number_ops_total; i++)
    {
      ni = variable_list.get_scalar_value(i);
      rhoi = variable_list.get_scalar_value(i + number_ops_total);
      ni = std::max(ni, constV(0.0));
      rhocalc += ni * ni * rhoi / sum_nsq;
    }

  for (unsigned int i = 0; i < number_ops_total; i++)
    {
      ni = variable_list.get_scalar_value(i);
      rhoi = variable_list.get_scalar_value(i + number_ops_total);
      ni = std::max(ni, constV(0.0));
      nix = variable_list.get_scalar_gradient(i);

      fnV = -ni + ni * ni * ni;
      for (unsigned int j = 0; j < number_ops_total; j++)
        {
          if (i != j)
            {
              nj = variable_list.get_scalar_value(j);
              nj = std::max(nj, constV(0.0));
              fnV += constV(2.0 * alpha) * ni * nj * nj;
            }
        }
      fnV *= m0;
      fsV = constV(fsbs) * 2.0 * ni * (rhoi - rhocalc) / sum_nsq;

      // Uncomment the next line to remove dependency on stored energy
      // fsV = constV(0.0);
      
      if (this->currentTime < userInputs.nucleation_start_time - 0.5 * userInputs.dtValue)
        {
          if ((i >= number_ops_ini) && (i <= number_ops_total - 1))
            {
              value_terms[i] = ni - gamma[i] * constV(userInputs.dtValue * MnV) * fnV;
              gradient_terms[i] = gamma[i] * constV(-userInputs.dtValue * KnV * MnV) * nix;
            }
          else
            {
              value_terms[i] = ni - constV(userInputs.dtValue * MnV) * fnV;
              gradient_terms[i] = constV(-userInputs.dtValue * KnV * MnV) * nix;
            }
        }
      else if (this->currentTime >=
                 userInputs.nucleation_start_time - 0.5 * userInputs.dtValue &&
               this->currentTime <= userInputs.nucleation_end_time)
        {
          if ((i >= number_ops_ini) && (i <= number_ops_total - 1))
            {
              value_terms[i] = ni + source_term[i];
              gradient_terms[i] = nix * constV(0.0);
            }
          else
            {
              value_terms[i] = ni;
              gradient_terms[i] = nix * constV(0.0);
            }
        }
      else if (this->currentTime <= userInputs.nucleation_end_time + globalSeedingTime)
        {
          if ((i >= number_ops_ini) && (i <= number_ops_total - 1))
            {
              value_terms[i] = ni - gamma[i] * constV(userInputs.dtValue * MnV) * fnV;
              gradient_terms[i] = gamma[i] * constV(-userInputs.dtValue * KnV * MnV) * nix;
            }
          else
            {
              value_terms[i] = ni - constV(userInputs.dtValue * MnV) * fnV;
              gradient_terms[i] = constV(-userInputs.dtValue * KnV * MnV) * nix;
            }
        }
      else
        {
          if ((i >= number_ops_ini) && (i <= number_ops_total - 1))
            {
              value_terms[i] = ni - gamma[i]*constV(userInputs.dtValue * MnV) * (fnV + fsV);
              gradient_terms[i] = gamma[i] * constV(-userInputs.dtValue * KnV * MnV) * nix;
            }
          else
            {
              value_terms[i] = ni - constV(userInputs.dtValue * MnV) * (fnV + fsV);
              gradient_terms[i] = constV(-userInputs.dtValue * KnV * MnV) * nix;
            }
        }
    }

  // Set the values of the RHS for output

  for (unsigned int i = 0; i < number_ops_total; i++)
    {
      variable_list.set_scalar_value_term_RHS(i, value_terms[i]);
      variable_list.set_scalar_gradient_term_RHS(i, gradient_terms[i]);
    }

  for (unsigned int i = 0; i < N_rho; i++)
    {
      rhoi = variable_list.get_scalar_value(i + number_ops_total);
      variable_list.set_scalar_value_term_RHS(i + number_ops_total, rhoi);
      variable_list.set_scalar_gradient_term_RHS(i + number_ops_total, nix * constV(0.0));
    }

  variable_list.set_scalar_value_term_RHS(number_ops_total + N_rho, rhocalc);

}

// seedNucleus
// a function particular to this app
template <int dim, int degree>
void
customPDE<dim, degree>::seedNucleus(
  const Point<dim, VectorizedArray<double>> &q_point_loc,
  std::vector<VectorizedArray<double>>      &source_term,
  std::vector<VectorizedArray<double>>      &gamma) const
{
  for (typename std::vector<nucleus<dim>>::const_iterator thisNucleus =
         this->nuclei.begin();
       thisNucleus != this->nuclei.end();
       ++thisNucleus)
    {
      if (this->currentTime < userInputs.nucleation_end_time + thisNucleus->seedingTime)
        {
          VectorizedArray<double> weighted_dist = this->weightedDistanceFromNucleusCenter(
            thisNucleus->center,
            userInputs.get_nucleus_freeze_semiaxes(thisNucleus->orderParameterIndex),
            q_point_loc,
            thisNucleus->orderParameterIndex);

          for (unsigned int i = 0; i < gamma[0].size(); i++)
            {
              if (weighted_dist[i] <= 1.0)
                {
                  gamma[thisNucleus->orderParameterIndex][i] = 0.0;

                  // Seed a nucleus if it was added to the list of nuclei this time step
                  if (thisNucleus->seedingTimestep == this->currentIncrement)
                    {
                      // Find the weighted distance to the outer edge of the nucleus and
                      // use it to calculate the order parameter source term
                      Point<dim, double> q_point_loc_element;
                      for (unsigned int j = 0; j < dim; j++)
                        {
                          q_point_loc_element(j) = q_point_loc(j)[i];
                        }
                      double r = this->weightedDistanceFromNucleusCenter(
                        thisNucleus->center,
                        userInputs.get_nucleus_semiaxes(thisNucleus->orderParameterIndex),
                        q_point_loc_element,
                        thisNucleus->orderParameterIndex);
                      double avg_semiaxis = 0.0;

                      for (unsigned int j = 0; j < dim; j++)
                        {
                          avg_semiaxis += thisNucleus->semiaxes[j];
                        }
                      avg_semiaxis /= dim;

                      source_term[thisNucleus->orderParameterIndex][i] =
                        0.5 *
                        (1.0 - std::tanh(avg_semiaxis * (r - 1.0) / interface_coeff));
                    }
                }
            }
        }
    }
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time
// independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are
// not explicit time-dependent equations. It takes "variable_list" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is
// given by "q_point_loc". The function outputs two terms to variable_list --
// one proportional to the test function and one proportional to the gradient of
// the test function. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // Number of order parameters
  unsigned int number_ops_ini = 23;

  // Number of empty OPs (for nucleation)
  unsigned int number_ops_nuc = 60;

  // Number of dislocation density fields
  unsigned int N_rho = 83;

  // Total number of order parameters
  unsigned int number_ops_total = number_ops_ini + number_ops_nuc;

  scalarvalueType sum_nsq = constV(0.0);
  scalarvalueType ni = constV(0.0);
  scalarvalueType nj = constV(0.0);

  // Calculation of double product and triple product
  scalarvalueType mult2Op_calc = constV(0.0);

  for (unsigned int i = 0; i < number_ops_ini; i++)
    {
      ni = variable_list.get_scalar_value(i);
      for (unsigned int j = 0; j < number_ops_ini; j++)
        {
          nj = variable_list.get_scalar_value(j);
          if (j > i)
            {
              mult2Op_calc += ni * nj;
            }
        }
    }

  variable_list.set_scalar_value_term_RHS(number_ops_total + N_rho + 1, mult2Op_calc);
}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. The
// (x,y,z) location of that quadrature point is given by "q_point_loc". The
// function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function -- for the
// left-hand-side of the equation. The index for each variable in this list
// corresponds to the index given at the top of this file. If there are multiple
// elliptic equations, conditional statements should be sed to ensure that the
// correct residual is being submitted. The index of the field being solved can
// be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{}
