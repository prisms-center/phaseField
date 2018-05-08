// Unit tests for the class "FloodFiller"
#include "../../include/FloodFiller.h"

template <int dim>
class InitialConditionFloodFill : public dealii::Function<dim>
{
public:
  InitialConditionFloodFill () : dealii::Function<dim>(1){
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const dealii::Point<dim> &p, const unsigned int component=0) const{
      /*
      dealii::Point<2> center(0.4,0.4);
      if (p.distance(center) < 0.2){
          return 1.0;
      }
      else{
          return 0.0;
      }
      */

      double val;
      if (p[1] < 0.6 and p[0] > p[1]+1.0e-10){
          val = 1.0;
      }
      else if (p[0] < 0.1 and p[1] > 0.9){
          val = 1.0;
      }
      else {
          val = 0.0;
      }

      std::cout << p[0] << "\t" << p[1] << "\t" << val << std::endl;
      return val;
  };

};

template <int dim,typename T>
  bool unitTest<dim,T>::test_FloodFiller(){

    char buffer[100];

	std::cout << "\nTesting 'FloodFiller'... " << std::endl;

    bool pass = false;

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

    vectorType *solution_field;
    solution_field = new vectorType;
    matrixFreeObject.initialize_dof_vector(*solution_field,  0); *solution_field = 0;

    // Set the field value
    //solution_field->print(std::cout);
    VectorTools::interpolate (dof_handler, InitialConditionFloodFill<dim>(), *solution_field);
    solution_field->print(std::cout);

    // Create a FloodFiller object
    QGaussLobatto<dim> quadrature2 (degree+1);

    FloodFiller<dim, degree> test_object(triangulation, fe, quadrature2);
    std::vector<GrainSet<dim>> grain_sets;
    GrainSet<dim> test_grain_set;
    //grain_sets.push_back(test_grain_set);
    test_object.calcGrainSets(fe, dof_handler, solution_field, 0.1, grain_sets);

    /*
    std::cout << "For grain: " << grain_sets[0].getGrainIndex() << std::endl;

    std::vector<std::vector<dealii::Point<dim>>> vertex_list = grain_sets[0].getVertexList();

    for (unsigned int c=0; c< vertex_list.size(); c++){
        for (unsigned int v=0; v< dealii::Utilities::fixed_power<dim>(2.0); v++){
            std::cout << vertex_list[c][v] << "\t";
        }
        std::cout << std::endl;
    }
    */

    /*
    // Subtests
    unsigned int subtest_index = 0;
    bool result;
    bool pass = true;

    // Subtest 1
    subtest_index++;
    result = false;
    if (test_object.getMaxIterations() == 123){
        result = true;
    }
    std::cout << "Subtest " << subtest_index << " result for 'getMaxIterations': " << result << std::endl;

    pass = pass && result;
    */

	sprintf (buffer, "Test result for 'FloodFiller': %u\n", pass);
	std::cout << buffer;

	return pass;
}
