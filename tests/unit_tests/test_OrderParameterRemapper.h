#include "../../include/OrderParameterRemapper.h"

template <int dim>
class InitialConditionOrderParameterRemapper : public dealii::Function<dim>
{
public:
  InitialConditionOrderParameterRemapper (unsigned int _index) : dealii::Function<dim>(1), index(_index){
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const dealii::Point<dim> &p, const unsigned int component=0) const{

      double val;

      if (index == 0){
          dealii::Point<dim> center0(0.125, 0.125);
          dealii::Point<dim> center1(1.0, 1.0);

          if (p.distance(center0) < std::sqrt(2.0)*0.15 or p.distance(center1) < std::sqrt(2.0)*0.15){
              val = 1.0;
          }
          else {
              val = 0.0;
          }
      }
      else {
          val = 0.0;
      }

      std::cout << index << "\t" << p[0] << "\t" << p[1] << "\t" << val << std::endl;
      return val;
  };

private:
    unsigned int index;

};

template <int dim,typename T>
  bool unitTest<dim,T>::test_OrderParameterRemapper(){

    char buffer[100];

	std::cout << "\nTesting 'OrderParameterRemapper'... " << std::endl;

    bool pass = false;

    // -------- Create the test mesh and solution vectors

    // Create the test mesh
    parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);

    const unsigned int degree = 1;

    FESystem<dim> fe(FE_Q<dim>(QGaussLobatto<1>(degree+1)),1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);
    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);


    typename MatrixFree<dim,double>::AdditionalData additional_data;
    #if (DEAL_II_VERSION_MAJOR < 9 && DEAL_II_VERSION_MINOR < 5)
        additional_data.mpi_communicator = MPI_COMM_WORLD;
    #endif
    additional_data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::partition_partition;
    additional_data.mapping_update_flags = (update_values | update_gradients | update_JxW_values | update_quadrature_points);
    QGaussLobatto<1> quadrature (degree+1);

    ConstraintMatrix constraints;
    constraints.clear();

    dealii::MatrixFree<dim,double> matrixFreeObject;
    matrixFreeObject.clear();
    matrixFreeObject.reinit (dof_handler, constraints, quadrature, additional_data);

    vectorType *solution_field_0;
    solution_field_0 = new vectorType;
    matrixFreeObject.initialize_dof_vector(*solution_field_0,  0); *solution_field_0 = 0;

    // Set the value for field 0
    VectorTools::interpolate (dof_handler, InitialConditionOrderParameterRemapper<dim>(0), *solution_field_0);

    solution_field_0->update_ghost_values();

    vectorType *solution_field_1;
    solution_field_1 = new vectorType;
    matrixFreeObject.initialize_dof_vector(*solution_field_1,  0); *solution_field_1 = 0;

    // Set the value for field 1
    VectorTools::interpolate (dof_handler, InitialConditionOrderParameterRemapper<dim>(1), *solution_field_1);

    solution_field_1->update_ghost_values();

    std::vector<vectorType*> solution_fields;
    solution_fields.push_back(solution_field_0);
    solution_fields.push_back(solution_field_1);

    // ---------- Create the simplified grain representations -----------

    QGaussLobatto<dim> quadrature2 (degree+1);

    FloodFiller<dim, degree> test_object(fe, quadrature2);
    std::vector<GrainSet<dim>> grain_sets_0;
    test_object.calcGrainSets(fe, dof_handler, solution_field_0, 0.1, 0, grain_sets_0);

    std::vector<GrainSet<dim>> grain_sets_1;
    test_object.calcGrainSets(fe, dof_handler, solution_field_1, 0.1, 1, grain_sets_1);

    std::vector<GrainSet<dim>> grain_sets = grain_sets_0;
    grain_sets.insert(grain_sets.end(), grain_sets_1.begin(), grain_sets_1.end());

    std::cout << "num grains: " << grain_sets.size() << " from 0: " << grain_sets_0.size() << " from 1: " << grain_sets_1.size() << std::endl;

    for (unsigned int g=0; g<grain_sets_0.size(); g++){
        std::cout << "g0: " << g << "\t" << grain_sets_0.at(g).getOrderParameterIndex() << std::endl;
    }

    for (unsigned int g=0; g<grain_sets_1.size(); g++){
        std::cout << "g1: " << g << "\t" << grain_sets_1.at(g).getOrderParameterIndex() << std::endl;
    }

    for (unsigned int g=0; g<grain_sets.size(); g++){
        std::cout << "g: " << g << "\t" << grain_sets.at(g).getOrderParameterIndex() << std::endl;
    }

    std::vector<SimplifiedGrainRepresentation<dim>> simplified_grain_representations;
    for (unsigned int g=0; g<grain_sets.size(); g++){
        SimplifiedGrainRepresentation<dim> simplified_grain_representation(grain_sets.at(g));
        simplified_grain_representations.push_back(simplified_grain_representation);
    }

    for (unsigned int g=0; g<simplified_grain_representations.size(); g++){
        std::cout << "Grain: " << g << "\t" << simplified_grain_representations.at(g).getCenter() << "\t" << simplified_grain_representations.at(g).getRadius() << "\t" << simplified_grain_representations.at(g).getGrainId() << "\t" << simplified_grain_representations.at(g).getOrderParameterId() << "\t" << simplified_grain_representations.at(g).getOldOrderParameterId() << "\n";
    }

    std::vector<unsigned int> order_parameter_id_list;
    order_parameter_id_list.push_back(0);
    order_parameter_id_list.push_back(1);

    SimplifiedGrainManipulator<dim> simplified_grain_manipulator;
    simplified_grain_manipulator.reassignGrains(simplified_grain_representations, 0.6, order_parameter_id_list);

    for (unsigned int g=0; g<simplified_grain_representations.size(); g++){
        std::cout << "Grain: " << g << "\t" << simplified_grain_representations.at(g).getCenter() << "\t" << simplified_grain_representations.at(g).getRadius() << "\t" << simplified_grain_representations.at(g).getGrainId() << "\t" << simplified_grain_representations.at(g).getOrderParameterId() << "\t" << simplified_grain_representations.at(g).getOldOrderParameterId() << "\n";
    }

    // ---------- The actual test run of OrderParameterRemapper -----------

    solution_field_1->print(std::cout);

    typename DoFHandler<dim>::active_cell_iterator di = dof_handler.begin_active();

    while (di != dof_handler.end())
    {
        if (di->is_locally_owned()){
            std::vector<types::global_dof_index> dof_indices(4,0);
            di->get_dof_indices(dof_indices);
            for (unsigned int i=0; i < dof_indices.size(); i++){
                    std::cout << "ob " << solution_field_0->local_element(dof_indices.at(i)) << " " << solution_field_1->local_element(dof_indices.at(i)) <<  std::endl;
                }
            }
        ++di;
    }

    OrderParameterRemapper<dim> order_parameter_remapper;
    order_parameter_remapper.remap(simplified_grain_representations, solution_fields, dof_handler, fe.dofs_per_cell);

    di = dof_handler.begin_active();
    while (di != dof_handler.end())
    {
        if (di->is_locally_owned()){
            std::vector<types::global_dof_index> dof_indices(4,0);
            di->get_dof_indices(dof_indices);
            for (unsigned int i=0; i < dof_indices.size(); i++){
                    std::cout << "oa " << solution_field_0->local_element(dof_indices.at(i)) << " " << solution_field_1->local_element(dof_indices.at(i)) <<  std::endl;
                }
            }
        ++di;
    }

    solution_field_1->print(std::cout);

	sprintf (buffer, "Test result for 'OrderParameterRemapper': %u\n", pass);
	std::cout << buffer;

	return pass;
}
