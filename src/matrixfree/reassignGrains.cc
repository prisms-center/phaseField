#include "../../include/matrixFreePDE.h"

#include "../../include/FloodFiller.h"
#include "../../include/SimplifiedGrainRepresentation.h"
#include "../../include/OrderParameterRemapper.h"

//vmult operation for LHS
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::reassignGrains () {

    if ( (currentIncrement%100==0) or currentIncrement == 0 ){

        //log time
        computing_timer.enter_section("matrixFreePDE: reassignGrains");

        pcout << "Reassigning grains..." << std::endl;

        // Create the simplified grain representations
        QGaussLobatto<dim> quadrature2 (degree+1);
        FloodFiller<dim, degree> test_object(*FESet.at(0), quadrature2);

        std::vector<GrainSet<dim>> grain_sets;
        std::vector<GrainSet<dim>> single_OP_grain_sets;

        for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
            // Eventually there should be a check to see if this is one of the shared order parameter grains
            std::vector<GrainSet<dim>> single_OP_grain_sets;
            test_object.calcGrainSets(*FESet.at(0), *dofHandlersSet_nonconst.at(0), solutionSet.at(fieldIndex), 1.0e-4, fieldIndex, single_OP_grain_sets);

            grain_sets.insert(grain_sets.end(), single_OP_grain_sets.begin(), single_OP_grain_sets.end());
            single_OP_grain_sets.clear();
        }

        // Set the grain indices to unique values (handling continuity from one time to another is handled elsewhere, outside of FloodFiller)
        for (unsigned int g=0; g<grain_sets.size(); g++){
            grain_sets.at(g).setGrainIndex(g);
        }

        std::vector<SimplifiedGrainRepresentation<dim>> simplified_grain_representations;
        for (unsigned int g=0; g<grain_sets.size(); g++){
            SimplifiedGrainRepresentation<dim> simplified_grain_representation(grain_sets.at(g));

            pcout << "Grain: " << simplified_grain_representation.getGrainId() << " " << simplified_grain_representation.getOrderParameterId() << " Center: " << simplified_grain_representation.getCenter()(0) << " " << simplified_grain_representation.getCenter()(1) << std::endl;

            simplified_grain_representations.push_back(simplified_grain_representation);
        }

        std::vector<unsigned int> order_parameter_id_list;
        for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
            // Eventually there should be a check to see if this is one of the shared order parameter grains
            order_parameter_id_list.push_back(fieldIndex);
        }

        double threshold = 10.0;

        SimplifiedGrainManipulator<dim> simplified_grain_manipulator;
        simplified_grain_manipulator.reassignGrains(simplified_grain_representations, threshold, order_parameter_id_list);

        OrderParameterRemapper<dim> order_parameter_remapper;
        order_parameter_remapper.remap(simplified_grain_representations, solutionSet, *dofHandlersSet_nonconst.at(0), FESet.at(0)->dofs_per_cell, threshold);

        pcout << "Reassigning grains completed." << std::endl << std::endl;

        //end log
        computing_timer.exit_section("matrixFreePDE: reassignGrains");
    }
}

#include "../../include/matrixFreePDE_template_instantiations.h"
