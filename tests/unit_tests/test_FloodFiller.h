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

      //std::cout << p[0] << "\t" << p[1] << "\t" << val << std::endl;
      return val;
  };

};

template <int dim>
void setExpectedVertexLists(std::vector<std::vector<dealii::Point<dim>>> & expected_vertex_list0, std::vector<std::vector<dealii::Point<dim>>> & expected_vertex_list1){
    // expected_vertex_list0
    {std::vector<dealii::Point<dim>> vertex_set(dealii::Utilities::fixed_power<dim>(2.0));
    {dealii::Point<dim> p(0.0,0.0); vertex_set[0] = p;}
    {dealii::Point<dim> p(0.25,0.0); vertex_set[1] = p;}
    {dealii::Point<dim> p(0.0,0.25); vertex_set[2] = p;}
    {dealii::Point<dim> p(0.25,0.25); vertex_set[3] = p;}
    expected_vertex_list0.push_back(vertex_set);}

    {std::vector<dealii::Point<dim>> vertex_set(dealii::Utilities::fixed_power<dim>(2.0));
    {dealii::Point<dim> p(0.25,0.0); vertex_set[0] = p;}
    {dealii::Point<dim> p(0.5,0.0); vertex_set[1] = p;}
    {dealii::Point<dim> p(0.25,0.25); vertex_set[2] = p;}
    {dealii::Point<dim> p(0.5,0.25); vertex_set[3] = p;}
    expected_vertex_list0.push_back(vertex_set);}

    {std::vector<dealii::Point<dim>> vertex_set(dealii::Utilities::fixed_power<dim>(2.0));
    {dealii::Point<dim> p(0.5,0.0); vertex_set[0] = p;}
    {dealii::Point<dim> p(0.75,0.0); vertex_set[1] = p;}
    {dealii::Point<dim> p(0.5,0.25); vertex_set[2] = p;}
    {dealii::Point<dim> p(0.75,0.25); vertex_set[3] = p;}
    expected_vertex_list0.push_back(vertex_set);}

    {std::vector<dealii::Point<dim>> vertex_set(dealii::Utilities::fixed_power<dim>(2.0));
    {dealii::Point<dim> p(0.75,0.0); vertex_set[0] = p;}
    {dealii::Point<dim> p(1.0,0.0); vertex_set[1] = p;}
    {dealii::Point<dim> p(0.75,0.25); vertex_set[2] = p;}
    {dealii::Point<dim> p(1.0,0.25); vertex_set[3] = p;}
    expected_vertex_list0.push_back(vertex_set);}

    {std::vector<dealii::Point<dim>> vertex_set(dealii::Utilities::fixed_power<dim>(2.0));
    {dealii::Point<dim> p(0.25,0.25); vertex_set[0] = p;}
    {dealii::Point<dim> p(0.5,0.25); vertex_set[1] = p;}
    {dealii::Point<dim> p(0.25,0.5); vertex_set[2] = p;}
    {dealii::Point<dim> p(0.5,0.5); vertex_set[3] = p;}
    expected_vertex_list0.push_back(vertex_set);}

    {std::vector<dealii::Point<dim>> vertex_set(dealii::Utilities::fixed_power<dim>(2.0));
    {dealii::Point<dim> p(0.5,0.25); vertex_set[0] = p;}
    {dealii::Point<dim> p(0.75,0.25); vertex_set[1] = p;}
    {dealii::Point<dim> p(0.5,0.5); vertex_set[2] = p;}
    {dealii::Point<dim> p(0.75,0.5); vertex_set[3] = p;}
    expected_vertex_list0.push_back(vertex_set);}

    {std::vector<dealii::Point<dim>> vertex_set(dealii::Utilities::fixed_power<dim>(2.0));
    {dealii::Point<dim> p(0.75,0.25); vertex_set[0] = p;}
    {dealii::Point<dim> p(1.0,0.25); vertex_set[1] = p;}
    {dealii::Point<dim> p(0.75,0.5); vertex_set[2] = p;}
    {dealii::Point<dim> p(1.0,0.5); vertex_set[3] = p;}
    expected_vertex_list0.push_back(vertex_set);}

    {std::vector<dealii::Point<dim>> vertex_set(dealii::Utilities::fixed_power<dim>(2.0));
    {dealii::Point<dim> p(0.5,0.5); vertex_set[0] = p;}
    {dealii::Point<dim> p(0.75,0.5); vertex_set[1] = p;}
    {dealii::Point<dim> p(0.5,0.75); vertex_set[2] = p;}
    {dealii::Point<dim> p(0.75,0.75); vertex_set[3] = p;}
    expected_vertex_list0.push_back(vertex_set);}

    {std::vector<dealii::Point<dim>> vertex_set(dealii::Utilities::fixed_power<dim>(2.0));
    {dealii::Point<dim> p(0.75,0.5); vertex_set[0] = p;}
    {dealii::Point<dim> p(1.0,0.5); vertex_set[1] = p;}
    {dealii::Point<dim> p(0.75,0.75); vertex_set[2] = p;}
    {dealii::Point<dim> p(1.0,0.75); vertex_set[3] = p;}
    expected_vertex_list0.push_back(vertex_set);}

    // expected_vertex_list1
    {std::vector<dealii::Point<dim>> vertex_set(dealii::Utilities::fixed_power<dim>(2.0));
    {dealii::Point<dim> p(0.0,0.75); vertex_set[0] = p;}
    {dealii::Point<dim> p(0.25,0.75); vertex_set[1] = p;}
    {dealii::Point<dim> p(0.0,1.0); vertex_set[2] = p;}
    {dealii::Point<dim> p(0.25,1.0); vertex_set[3] = p;}

    expected_vertex_list1.push_back(vertex_set);}
}

template <int dim>
bool compareUnsortedVectors(std::vector<std::vector<dealii::Point<dim>>> vec1, std::vector<std::vector<dealii::Point<dim>>> vec2, double tolerance){
    if (vec1.size() != vec2.size()){
        return false;
    }

    for (unsigned int i=0; i<vec1.size(); i++){
        bool element_found = false;
        for (unsigned int j=0; j<vec2.size(); j++){
            bool elements_equal = true;
            for (unsigned int v=0; v < dealii::Utilities::fixed_power<dim>(2.0); v++){
                if (vec1[i][v].distance(vec2[j][v]) > tolerance){
                    elements_equal = false;
                    break;
                }
            }
            element_found = elements_equal;
            if (element_found){
                break;
            }
        }
        if (!element_found){
            return false;
        }
    }

    return true;
}

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

    AffineConstraints<double> constraints;
    constraints.clear();

    dealii::MatrixFree<dim,double> matrixFreeObject;
    matrixFreeObject.clear();
      #if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR < 4)
   matrixFreeObject.reinit (dof_handler, constraints, quadrature, additional_data);
#else
   matrixFreeObject.reinit (MappingFE< dim, dim >(FE_Q<dim>(QGaussLobatto<1>(degree+1))),
       dof_handler, constraints, quadrature, additional_data);
#endif

    vectorType *solution_field;
    solution_field = new vectorType;
    matrixFreeObject.initialize_dof_vector(*solution_field,  0); *solution_field = 0;

    // Set the field value
    VectorTools::interpolate (dof_handler, InitialConditionFloodFill<dim>(), *solution_field);

    solution_field->update_ghost_values();

    //solution_field->print(std::cout);

    // Create a FloodFiller object
    QGaussLobatto<dim> quadrature2 (degree+1);

    FloodFiller<dim, degree> test_object(fe, quadrature2);
    std::vector<GrainSet<dim>> grain_sets;
    test_object.calcGrainSets(fe, dof_handler, solution_field, 0.1, 1.1, 0, grain_sets);

    std::vector<std::vector<dealii::Point<dim>>> expected_vertex_list0, expected_vertex_list1;

    setExpectedVertexLists(expected_vertex_list0, expected_vertex_list1);


    bool result = false;
    bool result0, result1;
    if (grain_sets.size() == 2){

            std::vector<std::vector<dealii::Point<dim>>> vertex_list0 = grain_sets[0].getVertexList();

            std::vector<std::vector<dealii::Point<dim>>> vertex_list1 = grain_sets[1].getVertexList();

            // Get the order in the canonical order
            if (vertex_list0.size() == 1){
                vertex_list0.swap(vertex_list1);
            }

            if (vertex_list0.size() == 9){
                result0 = compareUnsortedVectors(vertex_list0, expected_vertex_list0, 1.0e-3);
                std::cout << "Subtest result for grain 0: " << result0 << std::endl;
            }
            else {
                result0 = false;
            }

            if (vertex_list1.size() == 1){
                result1 = compareUnsortedVectors(vertex_list1, expected_vertex_list1, 1.0e-3);
                std::cout << "Subtest result for grain 1: " << result1 << std::endl;
            }
            else {
                result1 = false;
            }


            result = result0 and result1;
    }

    pass = result;

	snprintf(buffer, sizeof(buffer), "Test result for 'FloodFiller': %u\n", pass);
	std::cout << buffer;

	return pass;
}
