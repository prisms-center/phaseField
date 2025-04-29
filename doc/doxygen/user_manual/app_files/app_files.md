# Other App Files {#app_files}

This section details the files that define each application. The norm is that substantial modifications to these files consistute a new application (in contrast to simply making changes to the input file ('parameters.in'). Modifying these files may require some knowledge of C++.

A very brief description of the purpose of each of these files is below, an in-depth discussion can be found in the following subsections:
- equations.cc: Specifies the attributes of the model variables and the model's residual equations
- ICs\_and\_BCs.cc: Specifies the initial conditions for each variable (unless the initial condition is read from file) and non-uniform Dirichlet boundary conditions, if applicable.
- postprocess.cc: Specifies variable other than the primary model variables to output, as well as the expressions to derive these variables.
- nucleation.cc: Contains the function that determines the probability density of nucleation.
- customPDE.h: Contains the prototypes of the functions for the application, also contains declarations of the model constants given in the input file.
- main.cc: Main C++ function that controls the flow of the simulation. Identical for all the example applications and it is unlikely that users will need to modify it.

In all of these files, the user can access user inputs from parameters.in via the userInputs object (e.g. the domain size, the time step size). See the documentation entry for userInputParameters for a list of variable names inside ''userInputs''.

## equations.cc
The file ''equations.cc'' contains a list of the variables in the model equations and their attributes as well as the residuals for the model equations. The file contains four functions: loadVariableAttributes, explicitEquationRHS, nonExplicitEquationRHS, and equationLHS.

To modify the functions in this file, one needs to be familiar with the weak form of the governing equations. In PRISMS-PF, the governing equations are expressed in two terms. The first is the part of the integrand that is multiplied by the test function (marked by 'eq' with the subscript of the variable in the example below). The second is the part of the integrand that multiplied by the gradient of the test function (marked by 'eqx' with the subscript of the variable in the example below). For the coupled Cahn-Hilliard/Allen-Cahn system, the governing equations are

$$
\begin{align}
\int_{\Omega}   w  \eta^{n+1}  ~dV =\int_{\Omega}  w  \left( \eta^{n} - \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n) \right)+ \nabla w \cdot (- \Delta t M_{\eta}\kappa) \nabla \eta^{n} ~dV
\end{align}
$$

$$
\begin{align}
eq_{\eta} &= \eta^{n} - \Delta t M_{\eta}~ ((f_{\beta,c}^n-f_{\alpha,c}^n)H_{,\eta}^n)
\end{align}
$$

$$
\begin{align}
eqx_{\eta} &= (- \Delta t M_{\eta}\kappa) \nabla \eta^{n}
\end{align}
$$

and

$$
\begin{align}
]$\int_{\Omega}   w  c^{n+1}  ~dV = \int_{\Omega}   w c^{n} +  \nabla w   (-\Delta t M_{c})~ [~(f_{\alpha,cc}^n(1-H^{n+1})+f_{\beta,cc}^n H^{n+1}) \nabla c + ~((f_{\beta,c}^n-f_{\alpha,c}^n)H^{n+1}_{,\eta} \nabla \eta) ] ~dV
\end{align}
$$

$$
\begin{align}
eq_c &= c^{n}
\end{align}
$$

$$
\begin{align}
eqx_{c} &= (-\Delta t M_{c})~ [~(f_{\alpha,cc}^n(1-H^{n+1})+f_{\beta,cc}^n H^{n+1}) \nabla c + ~((f_{\beta,c}^n-f_{\alpha,c}^n)H^{n+1}_{,\eta} \nabla \eta)
\end{align}
$$

for the Allen-Cahn and Cahn-Hilliard equation, respectively. Each of the terms in the governing equation is marked with an underbrace. The terms multiplied by the test function are referred to as the value terms and the terms multiplied by the gradient of the test function are referred to as the gradient terms.

### loadVariableAttributes
Here is the loadVariableAttributes function for the coupled Allen-Cahn/Cahn-Hilliard example application:
```
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"c");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "c");
    set_dependencies_gradient_term_RHS(0, "n,grad(c)");

    // Variable 1
	set_variable_name				(1,"n");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(1, "c,n");
    set_dependencies_gradient_term_RHS(1, "grad(n)");
}
```

This function specifies the model variables and their attributes. In this case, the two model variables are the concentration, **c**, and the order parameter, **n**. Here, **c** is listed as the zeroth variable and **n** is listed as the first. For each variable, a series of attributes are set using a series of C++ function calls. The following table lists the functions and a description of the attributes they set:

| Function          | Options | Required | Default | Description |
| --------------|---------|----------|---------|----------------------------------------------------|
set_variable_name | [String] | no | var  | Sets the name of the variable. This name is used in 'parameters.in' as well as during output.
set_variable_type | SCALAR, VECTOR | no | SCALAR  | Sets whether the variable is a scalar or a vector.
set_variable_equation_type | EXPLICIT_TIME_DEPENDENT, AUXILIARY, TIME_INDEPENDENT | no | EXPLICIT_TIME_DEPENDENT  | Sets whether the governing equation for the variable is a time-dependent PDE (EXPLICIT_TIME_DEPENDENT), a time-independent PDE that does not require a linear solve (AUXILIARY) or a time independent PDE that does require a (non)linear solve (TIME_INDEPENDENT).
set_dependencies_value_term_RHS | String | yes | N/A| Sets which variables and their derivatives are needed to calculate the value term for the RHS. Variables are referenced by their names. First derivatives are referenced by ```grad``` and then the variable name in parentheses. Second derivatives are referenced by ```hess``` and then the variable name in parentheses.
set_dependencies_gradient_term_RHS | String | yes | N/A | Sets which variables and their derivatives are needed to calculate the gradient term for the RHS. Variables are referenced by their names. First derivatives are referenced by ```grad``` and then the variable name in parentheses. Second derivatives are referenced by ```hess``` and then the variable name in parentheses.
set_dependencies_value_term_LHS | String | no | [empty] | Sets which variables and their derivatives are needed to calculate the value term for the RHS. Variables are referenced by their names. First derivatives are referenced by ```grad``` and then the variable name in parentheses. Second derivatives are referenced by ```hess``` and then the variable name in parentheses. (Only needed for TIME_INDEPENDENT equations.)
set_dependencies_gradient_term_LHS | String | no | [empty] | Sets which variables and their derivatives are needed to calculate the gradient term for the RHS. Variables are referenced by their names. First derivatives are referenced by ```grad``` and then the variable name in parentheses. Second derivatives are referenced by ```hess``` and then the variable name in parentheses. (Only needed for TIME_INDEPENDENT equations.)
set_allowed_to_nucleate | Boolean | no | false | Sets whether the nucleation algorithms should be activated for this variable. (Only needed when nucleation is desired).
set_need_value_nucleation | Boolean | no | false | Sets whether the value of the variable is needed to calculate the nucleation probability in the 'nucleation.cc' file. (Only needed when nucleation is desired).

Some of these function calls are not present in the 'equations.cc' file for the coupledCahnHilliardAllenCahn application. Use of the LHS function calls can be found in the preciptiateEvolution application (among others) and use the nucleation function calls can be found the nucleationModel and nucleationModel_preferential applications.

### explicitEquationRHS
The explicitEquationRHS function is where the terms in the RHS of the governing equations for EXPLICIT_TIME_DEPENDENT equations are entered. The terms in the RHS of other equations are entered into the nonExplicitEquationRHS function. Here is the explicitEquationRHS function from the coupled Allen-Cahn/Cahn-Hilliard example application:
```
template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---

//c
scalarvalueType c = variable_list.get_scalar_value(0);
scalargradType cx = variable_list.get_scalar_gradient(0);

//n
scalarvalueType n = variable_list.get_scalar_value(1);
scalargradType nx = variable_list.get_scalar_gradient(1);

// --- Setting the expressions for the terms in the governing equations ---

// Free energy for each phase and their first and second derivatives
scalarvalueType fa = (-1.6704-4.776*c+5.1622*c*c-2.7375*c*c*c+1.3687*c*c*c*c);
scalarvalueType fac = (-4.776 + 10.3244*c - 8.2125*c*c + 5.4748*c*c*c);
scalarvalueType facc = (10.3244-16.425*c+16.4244*c*c);
scalarvalueType fb = (5.0*c*c-5.9746*c-1.5924);
scalarvalueType fbc = (10.0*c-5.9746);
scalarvalueType fbcc = constV(10.0);

// Interpolation function and its derivative
scalarvalueType h = (10.0*n*n*n-15.0*n*n*n*n+6.0*n*n*n*n*n);
scalarvalueType hn = (30.0*n*n-60.0*n*n*n+30.0*n*n*n*n);

// Residual equations
scalargradType mux = ( cx*((1.0-h)*facc+h*fbcc) + nx*((fbc-fac)*hn) );
scalarvalueType eq_c = c;
scalargradType eqx_c = (constV(-Mc*userInputs.dtValue)*mux);
scalarvalueType eq_n = (n-constV(userInputs.dtValue*Mn)*(fb-fa)*hn);
scalargradType eqx_n = (constV(-userInputs.dtValue*Kn*Mn)*nx);


// --- Submitting the terms for the governing equations ---

// Terms for the equation to evolve the concentration
variable_list.set_scalar_value_term_RHS(0,eq_c);
variable_list.set_scalar_gradient_term_RHS(0,eqx_c);

// Terms for the equation to evolve the order parameter
variable_list.set_scalar_value_term_RHS(1,eq_n);
variable_list.set_scalar_gradient_term_RHS(1,eqx_n);

}
```

In this function the equation terms at a particular quadrature point are calculated. The inputs to this function are a list of the model variable values and derivatives, variable_list and a point giving access to (x,y,z) coordinates, q_point_loc. The equation terms are added to variable_list as the output. The first few lines of the function set more convenient names for the variables and their derivatives. By convention, the value of the variable is denoted by the variable name (c for the concentration and n for the structural order parameter in this case), the list of first derivatives is denoted by the variable name followed by an ''x'', and second derivatives are denoted by the variable name followed by ''xx''. Each variable in variable can be accessed by the index it was given in loadVariableAttributes. The variable value and the derivatives can be accessed through the get_scalar_value,  get_scalar_gradient, and get_scalar_hessian object members for scalar variables and the  get_vector_value,  get_vector_gradient, and get_vector_hessian functions. The data type for the value of a scalar variable is scalarvalueType (a scalar), the data type for the first derivatives of a scalar variable is scalargradType (a vector with a length equal to the number of dimensions), and the data type for the second derivatives of a scalar variable is scalarhessType (a matrix with a size equal to the number of dimensions by the number of dimensions). For vector variables, the data types are vectorvalueType (a vector with length equal to the number of dimensions), vectorgradType (a matrix with a size equal to the number of dimensions by the number of dimensions), and vectorgradType (a rank-three tensor with a size in each direction equal to the number of dimensions).

After the nicknames for the field variables are set, the equation terms are calculated (including some intermediate variables, such as the free energies and the interpolation functions in the example above). These use the same six data types discussed in the preceding paragraph. The model variables given in 'parameters.in' can be used to define the residual terms.

Finally, the terms for the governing equations are submitted to variable_list, using the functions set_scalar_value_term_RHS, set_scalar_gradient_term_RHS, set_vector_value_term_RHS, and set_vector_gradient_term_RHS, once again referring the variables by their index.

The nickname step can be skipped, if desired, although the code is generally more readable with nicknames (and the performance is in some cases better using the nicknames). An example without nicknames for the fields can be found in the grainGrowth application, where loops over the ten variables are easier to construct when not declaring all of the variables at once. One word of caution, though: calling set_scalar_value_term_RHS and set_scalar_gradient_term_RHS overwrites the value and gradient of the variable, respectively (and same for the variants for vectors). Thus, either the equation terms need to be cached and then set after all the equation terms have been calculated (as is done in the grainGrowth application) or the values and gradients need to be cached before setting the equation terms (as is done with the standard nicknaming approach).

#### A note on types
The deal.II library uses a data structure called a VectorizedArray to store the variable values and their derivatives. This data structured in optimized for modern vectorized processors, giving a substantial speedup in some cases. However, this data structure can complicate things slightly. One complicating factor is that VectorizedArrays can't always be added, subtracted, multiplied, or divided with more standard data types like doubles. For this reason, you will see the ''constV(argument)'' function scattered throughout the code. This function turns a non-VectorizedArray into a VectorizedArray. To be safe, you can always encase non-VectorizedArrays with ''constV()'' when they share an operation with a VectorizedArray. A second complication is that not all of the standard mathematical operations are available for VectorizedArrays. The basic trigonometric functions are available, as are exponentials and square roots. However, hyperbolic tangents are not. If needed, they must be constructed from exponents (or by iterating through the VectorizedArray, see below). A third complication is that conditional statements involving VectorizedArrays are not allowed. To perform a conditional statement, you must iterate through the VectorizedArray. To do still, construct a for loop where the maximum index is [variable name].n_array_elements. For examples of this, refer to the postprocessing file for the grainGrowth app or the seedNucleus function in 'equations.cc' in the nucleationModel app. For more details on the deal.II implementation of VectorizedArrays (including a list of mathematical operations that are allowed), please visit [the relevant deal.II documentation page](https://www.dealii.org/8.4.0/doxygen/deal.II/classVectorizedArray.html).


### nonExplicitEquationRHS
The 'nonExplicitEquationRHS' function is where the terms in the RHS of the governing equations for AUXILIARY and TIME_INDEPENDENT equations are entered. The terms in the RHS of EXPLICIT_TIME_DEPENDENT equations are entered into the explicitEquationRHS function. The structure and use is otherwise identical to explicitEquationRHS. Examples of apps where this function is used include cahnHilliard (for an AUXILIARY equation) and MgNd_precipitate_single_Bppp (for AUXILIARY and TIME_INDEPENDENT equations).

### equationLHS
In the coupledCahnHilliardAllenCahn app, the equationLHS function is empty because it is only needed for TIME_INDEPENDENT PDEs (or more specifically, when a non-trivial matrix inversion needs to be performed). Here we go through the equationLHS function from the precipitateEvolution application, where the equation for mechanical equilibrium is TIME_INDEPENDENT.

From the formulation file in the precipitateEvolution application, the governing equation for the mechanical displacement is:

$$
\begin{align}
R(u) = \int_{\Omega}   \nabla w :  C(\eta_1, \eta_2, \eta_3) : \left( \epsilon - \epsilon^0(c,\eta_1, \eta_2, \eta_3)\right) ~dV = 0
\end{align}
$$


In PRISMS-PF, matrix inversion problems are always written as Newton's method iterations. For linear equations, like the one above, the solution is reached in a single Newton step. The reason for this approach is two-fold. First, it provides an identical user interface for linear and nonlinear problems. Second, it enables the efficient handling of constraints for when inhomogeneous Dirichlet boundary conditions are used.

To write the above equations in terms of a Newton iteration, the solution, $u$, can be written as the sum of an initial guess, $u_0$, and an update, $\Delta u$:

$$
\begin{align}
R(u) = R(u_0 + \Delta u) = R(u_0) +  \int_{\Omega} \left. \frac{\delta R(u)}{\delta u}\right|_{u=u_0} \Delta u ~dV = 0
\end{align}
$$

In this case, the equation is linear and the variation derivative is trivial:

$$
\begin{align}
R(u_0 + \Delta u) =  \int_{\Omega}   \nabla w :  C(\eta_1, \eta_2, \eta_3) : \left( \epsilon(u_0 + \Delta u) - \epsilon^0(c,\eta_1, \eta_2, \eta_3)\right) ~dV = \int_{\Omega}   \nabla w :  C(\eta_1, \eta_2, \eta_3) : \left( \epsilon(u_0) + \epsilon(\Delta u) - \epsilon^0(c,\eta_1, \eta_2, \eta_3)\right) ~dV = 0
\end{align}
$$

Rearranging yields:

$$
\begin{align}
\int_{\Omega} \nabla w : C : \nabla (\epsilon(\Delta u)) dV = -\int_{\Omega}   \nabla w : C :(\epsilon(u_0)-\epsilon^0) ~dV
\end{align}
$$

$$
\begin{align}
eqx_{u}^{LHS} &= C : \nabla (\epsilon(\Delta u))
\end{align}
$$

$$
\begin{align}
eqx_{u}^{RHS} &= C :(\epsilon(u_0)-\epsilon^0)
\end{align}
$$

The above values of $eqx_{u}^{LHS}$ and $eqx_{u}^{RHS}$ are used to define the residuals in the equations.h file. A similar process can be undertaken for other TIME_INDEPENDENT problems.

Here is the 'equationLHS' function from the precipitateEvolution example application:
```
template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---

//n1
scalarvalueType n1 = variable_list.get_scalar_value(1);

//n2
scalarvalueType n2 = variable_list.get_scalar_value(2);

//n3
scalarvalueType n3 = variable_list.get_scalar_value(3);

//u
vectorgradType Dux = variable_list.get_change_in_vector_gradient(4);


// --- Setting the expressions for the terms in the governing equations ---

vectorgradType eqx_Du;

// Interpolation functions

scalarvalueType h1V = (10.0*n1*n1*n1-15.0*n1*n1*n1*n1+6.0*n1*n1*n1*n1*n1);
scalarvalueType h2V = (10.0*n2*n2*n2-15.0*n2*n2*n2*n2+6.0*n2*n2*n2*n2*n2);
scalarvalueType h3V = (10.0*n3*n3*n3-15.0*n3*n3*n3*n3+6.0*n3*n3*n3*n3*n3);

// Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the dealii "symmetrize" function
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E;
E = symmetrize(Dux);

// Compute stress tensor (which is equal to the residual, Rux)
if (n_dependent_stiffness == true){
	dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ_combined;
	CIJ_combined = CIJ_Mg*(constV(1.0)-h1V-h2V-h3V) + CIJ_Beta*(h1V+h2V+h3V);

	computeStress<dim>(CIJ_combined, E, eqx_Du);
}
else{
	computeStress<dim>(CIJ_Mg, E, eqx_Du);
}

// --- Submitting the terms for the governing equations ---

variable_list.set_vector_gradient_term_LHS(4,eqx_Du);

}
```

Like the functions for the RHS, equationLHS takes variables_list and q_point_loc as inputs. However, because only one TIME_INDEPENDENT equation is solved at a time, the output are the terms for a single governing equation. As in explicitEquationRHS and nonExplicitEquationRHS, the model variable values and derivatives are given convenient names at the start of the file. The equationLHS function has a new option for getting variable values and derivatives, here one can access the change in the value of one of the variables or its derivatives. The function name is ```get_change_in_vector_gradient``` above. The other options are ```get_change_in_scalar_value```, ```get_change_in_scalar_gradient```, ```get_change_in_scalar_hessian```, ```get_change_in_vector_value```, and ```get_change_in_vector_hessian```.  Similar to the RHS functions, the value of eqx_Du is set and then used to set the appropriate residual for variables_list (in this case the vector residual term for the fourth variable).

If multiple TIME_INDEPENDENT equations need to be solved, you will need to use conditional statements to calculate the proper residual depending on the elliptic equation being solved. The index of the field being solved can be accessed using the ```this-$>$currentFieldIndex``` statement.

## ICs_and_BCs.cc
The file ''ICs_and_BCs.cc'' contains the initial conditions for the primary variables as well as the expressions for non-uniform Dirichlet boundary conditions, if present. The file contains two functions: setInitialCondition and setNonUniformDirichletBCs.

### setInitialCondition
Here's a look at the setInitialCondition function for the coupled Allen-Cahn/Cahn-Hilliard example application:
```
template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index
	  // The initial condition is two circles/spheres defined
	  // by a hyperbolic tangent function. The center of each circle/sphere is
	  // given by "center" and its radius is given by "rad".

      double dist;
	  scalar_IC = 0;

	  if (index == 0){
		  scalar_IC = matrix_concentration;
	  }

      dist = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++){
          dist += (p[dir]-center1[dir])*(p[dir]-center1[dir]);
      }
      dist = std::sqrt(dist);

      // Initial condition for the concentration field
      if (index == 0){
          scalar_IC += 0.5*(0.125)*(1.0-std::tanh((dist-radius1)/(1.0)));
      }
      else {
          scalar_IC += 0.5*(1.0-std::tanh((dist-radius1)/(1.0)));
      }

      dist = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++){
          dist += (p[dir]-center2[dir])*(p[dir]-center2[dir]);
      }
      dist = std::sqrt(dist);

      // Initial condition for the concentration field
      if (index == 0){
          scalar_IC += 0.5*(0.125)*(1.0-std::tanh((dist-radius2)/(1.0)));
      }
      else {
          scalar_IC += 0.5*(1.0-std::tanh((dist-radius2)/(1.0)));
      }


	  // --------------------------------------------------------------------------
}
```

In this function, the variable ```scalarIC``` should be set to the desired initial condition. The (x,y,z) coordinates can be accessed by the ```p``` variable, with indices zero through 2. For example, if you wanted to initialize a field to \f$5xy\f$, you would enter:

```
scalarIC = 5.0 * p[0] * p[1];
```

The initial conditions given above are a bit more complicated. They set the initial conditions for two circular (spherical in 3D) particles. The initial condition for different variables is set using conditional statements and the ```index``` variable. The ```index``` variable refers to the variable index from the top of the ''equations.h''. This example uses conditional statements to set the initial conditions differently for the concentration (```index=0```) and the structural order parameter  (```index=1```). The initial conditions for both fields are calculated using hyperbolic tangent functions that make use of the distance between the centers of the particles and the point ```p```.

The initial conditions for vector fields are set using ```vector_IC```. As an example, here is a section of the setInitialCondition function in the MgNd_precipitate_single_Bppp app:
```
else {
          for (unsigned int d=0; d<dim; d++){
              vector_IC(d) = 0.0;
          }
      }
```

The model constants defined in parameters.in can be used in the initial condition function (or boundary condition function), just as in the equation functions in equations.cc.

### setNonUniformDirichletBCs
This function is needed when one or more boundaries are set to NON_UNIFORM_DIRICHLET. These functions are very similar to the initial condition functions above, the user sets an expression for ```scalar_BC``` or ```vector_BC``` at point ```p```. Once again, the variable ```index``` is used to differentiate between variables based on their index from the top of the 'equations.cc' file. For Dirichlet boundary conditions that vary in time, the current time can be accessed via the variable ```time```.

Currently, the one of the only two applications that use a non-uniform Dirichlet boundary condition is CHiMaD_benchmark6a (the other is CHiMaD\_benchmark6b). Here is the setNonUniformDirichletBCs function from CHiMaD_benchmark6a:
```
template <int dim, int degree>
void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
{
    // --------------------------------------------------------------------------
    // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the boundary condition for each variable
    // according to its variable index. This function can be left blank if there
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).

    if (index == 2){
        if (direction == 1){
            double x=p[0];
            double y=p[1];
            scalar_BC=std::sin(y/7.0);
        }
    }

    // -------------------------------------------------------------------------

}
```

## postprocess.cc
Unlike the files discussed so far in this section, the postprocess.cc file is optional. This file allows users to specify fields other than the primary model variable to output. It also includes an option to integrate over a field variable (e.g. to check conservation of a field). In the example problems, this file is often used to calculate the total free energy of the system, which should be monotonically decreasing for most phase field problems. It also is convenient for calculating somewhat complicated expressions such as chemical potentials, stresses, and strains.

Note: Two restrictions are currently in place for the postprocessor. First, postprocessing variables can only be scalars, not vectors. Second, the postprocessor won't work if the primary model variable with index zero is a vector. (This is because postprocessor uses the mesh of the zeroth primary variable to calculate and output the postprocessed field.) These restrictions will be lifted in future versions of PRISMS-PF.

The file contains two functions: loadPostProcessorVariableAttributes and postProcessedFields.

### loadPostProcessorVariableAttributes
Similar to the equations.cc file, postprocess.cc begins with a function called loadPostProcessorVariableAttributes, which sets the attributes of the postprocessing variables. Here is the loadPostProcessorVariableAttributes function from the coupledCahnHilliardAllenCahn application, where the postprocessed field is ```f_tot```, the total free energy:
```
void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);

    set_dependencies_value_term_RHS(0, "c,n,grad(n)");
    set_dependencies_gradient_term_RHS(0, "");

    set_output_integral         	(0,true);

}
```

| Function          | Options | Required | Default | Description |
| --------------|---------|----------|---------|----------------------------------------------------|
set_variable_name | String | no | var  | Sets the name of the variable.
set_variable_type | SCALAR | no | SCALAR  | Sets whether the variable is a scalar or a vector. Only SCALAR is currently allowed.
set_dependencies_value_term_RHS | String | yes | | Sets which variables and their derivatives are needed to calculate the value term for the RHS. Variables are referenced by their names. First derivatives are referenced by ```grad``` and then the variable name in parentheses. Second derivatives are referenced by ```hess``` and then the variable name in parentheses.
set_dependencies_gradient_term_RHS | String | yes | | Sets which variables and their derivatives are needed to calculate the gradient term for the RHS. Variables are referenced by their names. First derivatives are referenced by ```grad``` and then the variable name in parentheses. Second derivatives are referenced by ```hess``` and then the variable name in parentheses.
set_output_integral | Boolean | no | false | Sets whether the integral of the variable should be calculated and written to a file named 'integratedFields.txt'.

### postProcessedFields
The second function in the postprocess.h file is similar to the residual functions in equations.h. The primary difference is that, while the primary model fields are read in from ```variable_list```, the residuals are submitted to ```pp_variable_list```. Otherwise, the structure is similar and expressions can be copied from one function to another.

Here is the postProcessedFields from the coupledCahnHilliardAllenCahn application:
```
template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---

//c
scalarvalueType c = variable_list.get_scalar_value(0);

//n
scalarvalueType n = variable_list.get_scalar_value(1);
scalargradType nx = variable_list.get_scalar_gradient(1);

// --- Setting the expressions for the terms in the postprocessing expressions ---

// Free energy for each phase and their first and second derivatives
scalarvalueType fa = (-1.6704-4.776*c+5.1622*c*c-2.7375*c*c*c+1.3687*c*c*c*c);
scalarvalueType fb = (5.0*c*c-5.9746*c-1.5924);

// Interpolation function and its derivative
scalarvalueType h = (10.0*n*n*n-15.0*n*n*n*n+6.0*n*n*n*n*n);

// The homogenous free energy
scalarvalueType f_chem = (constV(1.0)-h)*fa + h*fb;

// The gradient free energy
scalarvalueType f_grad = constV(0.5*Kn)*nx*nx;

// The total free energy
scalarvalueType f_tot;
f_tot = f_chem + f_grad;

// --- Submitting the terms for the postprocessing expressions ---

pp_variable_list.set_scalar_value_term_RHS(0, f_tot);

}
```

## nucleation.cc
The nucleation.cc file is also an optional file that is only needed when explicit nucleation is needed for the application. It contains the function getNucleationProbability, which calculates the probability of nucleation in a given volume element ```dV``` at a particular point in space ```p``` at time ```this-$>$currentTime```. Please refer to Note 6 [on the input file page ](https://prisms-center.github.io/phaseField/doxygen_files/input_file.html) for a description of the nucleation model in PRISMS-PF.

The values of the primary field variables can be accessed through the ```variable_value``` input, using the variable index from equations.cc. Only variables where set_need_value_nucleation was set to true in the loadVariableAttributes function in equations.cc can be accessed.

Here is getNucleationProbability from the nucleationModel application:
```
template <int dim, int degree>
double customPDE<dim,degree>::getNucleationProbability(variableValueContainer variable_value, double dV, dealii::Point<dim> p, unsigned int variable_index) const
{
	//Supersaturation factor
    double ssf;
    if (dim ==2) ssf=variable_value(0)-calmin;
    if (dim ==3) ssf=(variable_value(0)-calmin)*(variable_value(0)-calmin);
	// Calculate the nucleation rate
	double J=k1*exp(-k2/(std::max(ssf,1.0e-6)))*exp(-tau/(this->currentTime));
	double retProb=1.0-exp(-J*userInputs.dtValue*((double)userInputs.steps_between_nucleation_attempts)*dV);
    return retProb;
}
```

## customPDE.h
The final file in each application folder that users often need to change is customPDE.h. This file contains the declarations for all of the functions and variables specific to the application. In C++ terminology, it is the declaration for the customPDE class, which is a subclass of MatrixFreePDE. For most users, the relevant section is labeled as ''Model constants specific to this subclass'', which is where the model constants from parameters.in are extracted from the ```userInputs``` object. There is a separate function to extract constants of each type. These are:

- get\_model\_constant\_double
- get\_model\_constant\_int
- get\_model\_constant\_bool
- get\_model\_constant\_rank\_1\_tensor
- get\_model\_constant\_rank\_2\_tensor
- get\_model\_constant\_elasticity\_tensor

Each of these take the variable name from parameters.in as an input and output the appropriate type.

Here is customPDE.h from the preciptiateEvolution application, where all of the different types of constants are used:
```
#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
	customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {};
    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);

    private:
    #include "../../include/typeDefs.h"

    const userInputParameters<dim> userInputs;

    // Function to set the RHS of the governing equations for explicit time dependent equations (in equations.h)
    void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the RHS of the governing equations for all other equations (in equations.h)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the LHS of the governing equations (in equations.h)
    void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set postprocessing expressions (in postprocess.h)
    #ifdef POSTPROCESS_FILE_EXISTS
    void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                    variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
                    const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
    #endif

    // Function to set the nucleation probability (in nucleation.h)
    #ifdef NUCLEATION_FILE_EXISTS
    double getNucleationProbability(variableValueContainer variable_value, double dV) const;
    #endif

	// ================================================================
	// Methods specific to this subclass
	// ================================================================


	// ================================================================
	// Model constants specific to this subclass
	// ================================================================

	double McV = userInputs.get_model_constant_double("McV");
	double Mn1V = userInputs.get_model_constant_double("Mn1V");
	double Mn2V = userInputs.get_model_constant_double("Mn2V");
	double Mn3V = userInputs.get_model_constant_double("Mn3V");
	dealii::Tensor<2,dim> Kn1 = userInputs.get_model_constant_rank_2_tensor("Kn1");
	dealii::Tensor<2,dim> Kn2 = userInputs.get_model_constant_rank_2_tensor("Kn2");
	dealii::Tensor<2,dim> Kn3 = userInputs.get_model_constant_rank_2_tensor("Kn3");
	bool n_dependent_stiffness = userInputs.get_model_constant_bool("n_dependent_stiffness");
	dealii::Tensor<2,dim> sfts_linear1 = userInputs.get_model_constant_rank_2_tensor("sfts_linear1");
	dealii::Tensor<2,dim> sfts_const1 = userInputs.get_model_constant_rank_2_tensor("sfts_const1");
	dealii::Tensor<2,dim> sfts_linear2 = userInputs.get_model_constant_rank_2_tensor("sfts_linear2");
	dealii::Tensor<2,dim> sfts_const2 = userInputs.get_model_constant_rank_2_tensor("sfts_const2");
	dealii::Tensor<2,dim> sfts_linear3 = userInputs.get_model_constant_rank_2_tensor("sfts_linear3");
	dealii::Tensor<2,dim> sfts_const3 = userInputs.get_model_constant_rank_2_tensor("sfts_const3");
	double A4 = userInputs.get_model_constant_double("A4");
	double A3 = userInputs.get_model_constant_double("A3");
	double A2 = userInputs.get_model_constant_double("A2");
	double A1 = userInputs.get_model_constant_double("A1");
	double A0 = userInputs.get_model_constant_double("A0");
	double B2 = userInputs.get_model_constant_double("B2");
	double B1 = userInputs.get_model_constant_double("B1");
	double B0 = userInputs.get_model_constant_double("B0");

	const static unsigned int CIJ_tensor_size =2*dim-1+dim/3;
	dealii::Tensor<2,CIJ_tensor_size> CIJ_Mg = userInputs.get_model_constant_elasticity_tensor("CIJ_Mg");
	dealii::Tensor<2,CIJ_tensor_size> CIJ_Beta = userInputs.get_model_constant_elasticity_tensor("CIJ_Beta");



	// ================================================================

};
```

This file can also be used to declare new member functions for the application. One example of this is the seedNucleus function in the nucleationModel application.

Furthermore, for advanced users, the customPDE class can be used to override MatrixFreePDE functions from the core PRISMS-PF library. One example of this is in the CHiMaD_benchmark6b application, where the makeTriangulation function is overridden to create a non-square mesh.

## main.cc
The final C++ file in the application directory is the main.cc file. This file controls the overall flow of the code and is unlikely to be modified by most users. For all of the example applications, main.cc is identical. One situation where a user may want to modify 'main.cc' is if they wanted to run several simulations with different parameter sets for one execution of the code.
