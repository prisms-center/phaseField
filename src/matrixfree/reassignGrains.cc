#include "../../include/matrixFreePDE.h"

#include "../../include/FloodFiller.h"
#include "../../include/SimplifiedGrainRepresentation.h"
#include "../../include/OrderParameterRemapper.h"

//vmult operation for LHS
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::reassignGrains () {

    //log time
    computing_timer.enter_section("matrixFreePDE: reassignGrains");

    pcout << "Reassigning grains..." << std::endl;

    // Get the index of the first scalar field (used to get the FE object and DOFHandler)
    unsigned int scalar_field_index = 0;
    for (unsigned int var=0; var<userInputs.number_of_variables; var++){
        if (userInputs.var_type == SCALAR){
            scalar_field_index = var;
            break;
        }
    }

    // Create the simplified grain representations
    QGaussLobatto<dim> quadrature2 (degree+1);
    FloodFiller<dim, degree> test_object(*FESet.at(scalar_field_index), quadrature2);

    std::vector<GrainSet<dim>> grain_sets;
    std::vector<GrainSet<dim>> single_OP_grain_sets;

    unsigned int op_list_index = 0;
    for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
        if (fieldIndex == userInputs.variables_for_remapping.at(op_list_index)){
            op_list_index++;

            // Eventually there should be a check to see if this is one of the shared order parameter variables
            std::vector<GrainSet<dim>> single_OP_grain_sets;
            test_object.calcGrainSets(*FESet.at(scalar_field_index), *dofHandlersSet_nonconst.at(scalar_field_index), solutionSet.at(fieldIndex), order_parameter_threshold, fieldIndex, single_OP_grain_sets);

            grain_sets.insert(grain_sets.end(), single_OP_grain_sets.begin(), single_OP_grain_sets.end());
            single_OP_grain_sets.clear();
        }
    }

    // Set the grain indices to unique values
    for (unsigned int g=0; g<grain_sets.size(); g++){
        grain_sets.at(g).setGrainIndex(g);
    }

    std::vector<SimplifiedGrainRepresentation<dim>> old_grain_representations = simplified_grain_representations;
    simplified_grain_representations.clear();
    for (unsigned int g=0; g<grain_sets.size(); g++){
        SimplifiedGrainRepresentation<dim> simplified_grain_representation(grain_sets.at(g));

        pcout << "Grain: " << simplified_grain_representation.getGrainId() << " " << simplified_grain_representation.getOrderParameterId() << " Center: " << simplified_grain_representation.getCenter()(0) << " " << simplified_grain_representation.getCenter()(1) << std::endl;

        simplified_grain_representations.push_back(simplified_grain_representation);
    }

    SimplifiedGrainManipulator<dim> simplified_grain_manipulator;

    if (currentIncrement > 0){
        simplified_grain_manipulator.transferGrainIds(old_grain_representations, simplified_grain_representations);
    }

    simplified_grain_manipulator.reassignGrains(simplified_grain_representations, userInputs.buffer_between_grains, variables_for_remapping);

    OrderParameterRemapper<dim> order_parameter_remapper;
    order_parameter_remapper.remap(simplified_grain_representations, solutionSet, *dofHandlersSet_nonconst.at(0), FESet.at(0)->dofs_per_cell, userInputs.buffer_between_grains);

    pcout << "Reassigning grains completed." << std::endl << std::endl;

    //end log
    computing_timer.exit_section("matrixFreePDE: reassignGrains");


}

#include "../../include/matrixFreePDE_template_instantiations.h"
