/*
 * customPDE.h
 *
 *  Created on: Feb 24, 2017
 *      Author: stephendewitt
 */

#ifndef APPLICATIONS_ALLENCAHN_CUSTOMPDE_H_
#define APPLICATIONS_ALLENCAHN_CUSTOMPDE_H_

#include "../../include/matrixFreePDE.h"
#include "../../include/parallelNucleationList.h"
#include "../../include/nucleus.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:

	// Constructor, which calls the MatrixFreePDE constructor
	customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) {};

private:
	#include "../../include/typeDefs.h"

	// Pure virtual method in MatrixFreePDE
	void residualRHS(const std::vector<modelVariable<dim> > & modelVarList,
			  	  	 std::vector<modelResidual<dim> > & modelResidualsList,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Pure virtual method in MatrixFreePDE
	void residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
	  		  	  	 modelResidual<dim> & modelRes,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Pure virtual method in MatrixFreePDE
	void energyDensity(const std::vector<modelVariable<dim> > & modelVarList, const dealii::VectorizedArray<double> & JxW_value,
			  	  	 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc);

	// Virtual method in MatrixFreePDE we choose to override
	void getNucleiList ();

	// ================================================================
	// Model constants
	// ================================================================

	double McV = boost::get<double>(this->userInputs.model_constants[0]);
	double MnV = boost::get<double>(this->userInputs.model_constants[1]);
	double KnV = boost::get<double>(this->userInputs.model_constants[2]);
	double c_avg = boost::get<double>(this->userInputs.model_constants[3]);
	double W_barrier = boost::get<double>(this->userInputs.model_constants[4]);
	double A0 = boost::get<double>(this->userInputs.model_constants[5]);
	double A2 = boost::get<double>(this->userInputs.model_constants[6]);
	double calmin = boost::get<double>(this->userInputs.model_constants[7]);
	double B0 = boost::get<double>(this->userInputs.model_constants[8]);
	double B2 = boost::get<double>(this->userInputs.model_constants[9]);
	double cbtmin = boost::get<double>(this->userInputs.model_constants[10]);

	double interface_coeff=std::sqrt(2.0*KnV/W_barrier);

	// ================================================================

	// ----------------------------------------------------------------
	// Nucleation methods specific to this subclass
	// ----------------------------------------------------------------

	// Function to determine where new nuclei are seeded, varies between applications unless generalized, accesses MatrixFreePDE members
	void getLocalNucleiList(std::vector<nucleus<dim> > &newnuclei) const;

	// Contains nucleation probability that varies between applications, no MatrixFreePDE member access
	double nucProb(double cValue, double dV, double ct) const;

    //Contains safety check to ensure that no new nuclei from complete prospective list overlap with existing precipitates
    void safetyCheckNewNuclei(std::vector<nucleus<dim> > newnuclei, std::vector<unsigned int> &conflict_inds);

	// Function to refine the mesh near the new nuclei, generic, accesses and modifies MatrixFreePDE members
	void refineMeshNearNuclei(std::vector<nucleus<dim> > newnuclei);

	// Vector of all the nuclei seeded in the problem
	std::vector<nucleus<dim> > nuclei;
};

#endif /* APPLICATIONS_ALLENCAHN_CUSTOMPDE_H_ */
