#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
    // Constructor
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


    void adaptiveRefine(unsigned int _currentIncrement);

	// ================================================================
	// Methods specific to this subclass
	// ================================================================

	// ================================================================
	// Model constants specific to this subclass
	// ================================================================

	double D = userInputs.get_model_constant_double("D");
	double W0 = userInputs.get_model_constant_double("W0");
	double delta = userInputs.get_model_constant_double("delta");
	double epsilonM = userInputs.get_model_constant_double("epsilonM");
	double theta0 = userInputs.get_model_constant_double("theta0");
	double mult = userInputs.get_model_constant_double("mult");

    // The coupling constant, determined from solvability theory
    double lambda = (D/0.6267/W0/W0);

	// ================================================================

};


#include <deal.II/distributed/grid_refinement.h>

//default implementation of adaptive mesh refinement
template <int dim, int degree>
void customPDE<dim,degree>::adaptiveRefine(unsigned int currentIncrement){
if (this->userInputs.h_adaptivity == true){
	if ( (currentIncrement == 0) ){
		this->computing_timer.enter_section("matrixFreePDE: AMR");
		unsigned int numDoF_preremesh = this->totalDOFs;
		for (unsigned int remesh_index=0; remesh_index < (this->userInputs.max_refinement_level-this->userInputs.min_refinement_level); remesh_index++){

			this->adaptiveRefineCriterion();
			this->refineGrid();
			this->reinit();

			// If the mesh hasn't changed from the previous cycle, stop remeshing
			if (this->totalDOFs == numDoF_preremesh) break;
			numDoF_preremesh = this->totalDOFs;
		}
		this->computing_timer.exit_section("matrixFreePDE: AMR");
	}
	else {

        double remesh_coeff = 1;
        if (currentIncrement > 40000){
            remesh_coeff = 2.5;
        }
        else if (currentIncrement > 20000){
            remesh_coeff = 1.5;
        }

        if ( (currentIncrement%(int(remesh_coeff * this->userInputs.skip_remeshing_steps))==0) ){

    		this->computing_timer.enter_section("matrixFreePDE: AMR");

    		// Apply constraints before remeshing
    		for(unsigned int fieldIndex=0; fieldIndex<this->fields.size(); fieldIndex++){
    			this->constraintsDirichletSet[fieldIndex]->distribute(*this->solutionSet[fieldIndex]);
    			this->constraintsOtherSet[fieldIndex]->distribute(*this->solutionSet[fieldIndex]);
    			this->solutionSet[fieldIndex]->update_ghost_values();
    		}
    		this->adaptiveRefineCriterion();
    		this->refineGrid();
    		this->reinit();
    		this->computing_timer.exit_section("matrixFreePDE: AMR");
        }
	}
}
}
