#include "../../include/SimplifiedGrainRepresentation.h"

template <int dim, typename T>
bool
unitTest<dim, T>::test_SimplifiedGrainManipulator_transferGrainIds()
{
  char buffer[100];

  std::cout << "\nTesting 'SimplifiedGrainManipulator::transferGrainIds'... "
            << std::endl;

  bool pass = true;

  // Subtest 1 (two grains both needing reassignment)
  {
    // Create grain 0_old with center (0.5, 0.5), radius 0.5*sqrt(2), and order
    // parameter 0
    GrainSet<dim> test_grain_set_0;
    {
      test_grain_set_0.setGrainIndex(0);
      test_grain_set_0.setOrderParameterIndex(0);

      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(0.0, 0.0);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(1.0, 0.0);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(0.0, 1.0);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(1.0, 1.0);
        vertex_set[3] = p;
      }

      test_grain_set_0.addVertexList(vertex_set);
    }

    SimplifiedGrainRepresentation<dim> simplified_grain_representation_0(
      test_grain_set_0);

    // Create grain 1_old with center (4.5, 0.5), radius 0.5*sqrt(2), and order
    // parameter 1
    GrainSet<dim> test_grain_set_1;
    {
      test_grain_set_1.setGrainIndex(1);
      test_grain_set_1.setOrderParameterIndex(1);

      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(4.0, 0.0);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(5.0, 0.0);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(4.0, 1.0);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(5.0, 1.0);
        vertex_set[3] = p;
      }

      test_grain_set_1.addVertexList(vertex_set);
    }

    SimplifiedGrainRepresentation<dim> simplified_grain_representation_1(
      test_grain_set_1);

    // Create grain 0_new with center (0.6, 0.6), radius 0.6*sqrt(2), and order
    // parameter 0
    GrainSet<dim> test_grain_set_0_new;
    {
      test_grain_set_0_new.setGrainIndex(5);
      test_grain_set_0_new.setOrderParameterIndex(0);

      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(0.0, 0.0);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(1.2, 0.0);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(0.0, 1.2);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(1.2, 1.2);
        vertex_set[3] = p;
      }

      test_grain_set_0_new.addVertexList(vertex_set);
    }

    SimplifiedGrainRepresentation<dim> simplified_grain_representation_0_new(
      test_grain_set_0_new);

    // Create grain 1_new with center (4.6, 0.6), radius 0.6*sqrt(2), and order
    // parameter 1
    GrainSet<dim> test_grain_set_1_new;
    {
      test_grain_set_1.setGrainIndex(3);
      test_grain_set_1.setOrderParameterIndex(1);

      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(4.0, 0.0);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(5.2, 0.0);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(4.0, 1.2);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(5.2, 1.2);
        vertex_set[3] = p;
      }

      test_grain_set_1_new.addVertexList(vertex_set);
    }

    SimplifiedGrainRepresentation<dim> simplified_grain_representation_1_new(
      test_grain_set_1_new);

    // Now build of the vectors of these grain representations
    std::vector<SimplifiedGrainRepresentation<dim>> old_grain_representations,
      new_grain_representations;
    old_grain_representations.push_back(simplified_grain_representation_0);
    old_grain_representations.push_back(simplified_grain_representation_1);
    new_grain_representations.push_back(simplified_grain_representation_0_new);
    new_grain_representations.push_back(simplified_grain_representation_1_new);

    // Now run the actual test
    std::cout << "Grain ids before transfer: "
              << new_grain_representations.at(0).getGrainId() << " "
              << new_grain_representations.at(1).getGrainId() << std::endl;

    SimplifiedGrainManipulator<dim> simplified_grain_manipulator;
    simplified_grain_manipulator.transferGrainIds(old_grain_representations,
                                                  new_grain_representations);

    std::cout << "Grain ids after transfer: "
              << new_grain_representations.at(0).getGrainId() << " "
              << new_grain_representations.at(1).getGrainId() << std::endl;

    bool result = false;
    if (new_grain_representations.at(0).getGrainId() == 0 and
        new_grain_representations.at(1).getGrainId() == 1)
      {
        result = true;
      }
    pass = pass & result;
  }

  snprintf(buffer,
           sizeof(buffer),
           "Test result for 'SimplifiedGrainManipulator::transferGrainIds': %u\n",
           pass);
  std::cout << buffer;

  return pass;
}

template <int dim, typename T>
bool
unitTest<dim, T>::test_SimplifiedGrainManipulator_reassignGrains()
{
  char buffer[100];

  std::cout << "\nTesting 'SimplifiedGrainManipulator::reassignGrains'... " << std::endl;

  bool pass = true;

  // Subtest 1 (three grains, two order parameters, one grain needs
  // reassignment)
  {
    // Create grain 0 with center (0.5, 0.5), radius 0.5*sqrt(2), and order
    // parameter 0
    GrainSet<dim> test_grain_set_0;
    {
      test_grain_set_0.setGrainIndex(0);
      test_grain_set_0.setOrderParameterIndex(0);

      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(0.0, 0.0);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(1.0, 0.0);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(0.0, 1.0);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(1.0, 1.0);
        vertex_set[3] = p;
      }

      test_grain_set_0.addVertexList(vertex_set);
    }

    SimplifiedGrainRepresentation<dim> simplified_grain_representation_0(
      test_grain_set_0);

    // Create grain 1 with center (4.5, 0.5), radius 0.5*sqrt(2), and order
    // parameter 1
    GrainSet<dim> test_grain_set_1;
    {
      test_grain_set_1.setGrainIndex(1);
      test_grain_set_1.setOrderParameterIndex(1);

      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(4.0, 0.0);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(5.0, 0.0);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(4.0, 1.0);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(5.0, 1.0);
        vertex_set[3] = p;
      }

      test_grain_set_1.addVertexList(vertex_set);
    }

    SimplifiedGrainRepresentation<dim> simplified_grain_representation_1(
      test_grain_set_1);

    // Create grain 2 with center (0.6, 0.6), radius 0.6*sqrt(2), and order
    // parameter 0
    GrainSet<dim> test_grain_set_2;
    {
      test_grain_set_2.setGrainIndex(2);
      test_grain_set_2.setOrderParameterIndex(0);

      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(0.0, 0.0);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(1.2, 0.0);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(0.0, 1.2);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(1.2, 1.2);
        vertex_set[3] = p;
      }

      test_grain_set_2.addVertexList(vertex_set);
    }

    SimplifiedGrainRepresentation<dim> simplified_grain_representation_2(
      test_grain_set_2);

    // Now build of the vectors of these grain representations
    std::vector<SimplifiedGrainRepresentation<dim>> grain_representations;
    grain_representations.push_back(simplified_grain_representation_0);
    grain_representations.push_back(simplified_grain_representation_1);
    grain_representations.push_back(simplified_grain_representation_2);

    std::vector<unsigned int> order_parameter_id_list;
    order_parameter_id_list.push_back(0);
    order_parameter_id_list.push_back(1);

    // Now run the actual test
    /*
    for (unsigned int g=0; g < grain_representations.size(); g++){
        std::cout << "Grain ops before reassignment: " <<
    grain_representations.at(g).getOrderParameterId() << std::endl;
    }
    */
    SimplifiedGrainManipulator<dim> simplified_grain_manipulator;
    simplified_grain_manipulator.reassignGrains(grain_representations,
                                                0.5,
                                                order_parameter_id_list);
    /*
    for (unsigned int g=0; g < grain_representations.size(); g++){
        std::cout << "Grain ops after reassignment: " <<
    grain_representations.at(g).getOrderParameterId() << std::endl;
    }
    */
    bool result = false;
    if (grain_representations.at(0).getOrderParameterId() == 1 and
        grain_representations.at(1).getOrderParameterId() == 1 and
        grain_representations.at(2).getOrderParameterId() == 0)
      {
        result = true;
      }

    snprintf(buffer,
             sizeof(buffer),
             "Subtest 1 result for 'SimplifiedGrainManipulator::reassignGrains': %u\n",
             result);
    std::cout << buffer;

    pass = pass & result;
  }

  // Subtest 2 (three grains, two order parameters, no grains need reassignment)
  {
    // Create grain 0 with center (0.5, 0.5), radius 0.5*sqrt(2), and order
    // parameter 0
    GrainSet<dim> test_grain_set_0;
    {
      test_grain_set_0.setGrainIndex(0);
      test_grain_set_0.setOrderParameterIndex(0);

      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(0.0, 0.0);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(1.0, 0.0);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(0.0, 1.0);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(1.0, 1.0);
        vertex_set[3] = p;
      }

      test_grain_set_0.addVertexList(vertex_set);
    }

    SimplifiedGrainRepresentation<dim> simplified_grain_representation_0(
      test_grain_set_0);

    // Create grain 1 with center (4.5, 0.5), radius 0.5*sqrt(2), and order
    // parameter 1
    GrainSet<dim> test_grain_set_1;
    {
      test_grain_set_1.setGrainIndex(1);
      test_grain_set_1.setOrderParameterIndex(1);

      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(4.0, 0.0);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(5.0, 0.0);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(4.0, 1.0);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(5.0, 1.0);
        vertex_set[3] = p;
      }

      test_grain_set_1.addVertexList(vertex_set);
    }

    SimplifiedGrainRepresentation<dim> simplified_grain_representation_1(
      test_grain_set_1);

    // Create grain 2 with center (2.5, 2.5), radius 0.5*sqrt(2), and order
    // parameter 0
    GrainSet<dim> test_grain_set_2;
    {
      test_grain_set_2.setGrainIndex(2);
      test_grain_set_2.setOrderParameterIndex(0);

      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(2.0, 2.0);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(3.0, 2.0);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(2.0, 3.0);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(3.0, 3.0);
        vertex_set[3] = p;
      }

      test_grain_set_2.addVertexList(vertex_set);
    }

    SimplifiedGrainRepresentation<dim> simplified_grain_representation_2(
      test_grain_set_2);

    // Now build of the vectors of these grain representations
    std::vector<SimplifiedGrainRepresentation<dim>> grain_representations;
    grain_representations.push_back(simplified_grain_representation_0);
    grain_representations.push_back(simplified_grain_representation_1);
    grain_representations.push_back(simplified_grain_representation_2);

    std::vector<unsigned int> order_parameter_id_list;
    order_parameter_id_list.push_back(0);
    order_parameter_id_list.push_back(1);

    // Now run the actual test
    /*
    for (unsigned int g=0; g < grain_representations.size(); g++){
        std::cout << "Grain ops before reassignment: " <<
    grain_representations.at(g).getOrderParameterId() << std::endl;
    }
    */

    SimplifiedGrainManipulator<dim> simplified_grain_manipulator;
    simplified_grain_manipulator.reassignGrains(grain_representations,
                                                0.5,
                                                order_parameter_id_list);

    /*
    for (unsigned int g=0; g < grain_representations.size(); g++){
        std::cout << "Grain ops after reassignment: " <<
    grain_representations.at(g).getOrderParameterId() << std::endl;
    }
    */

    bool result = false;
    if (grain_representations.at(0).getOrderParameterId() == 0 and
        grain_representations.at(1).getOrderParameterId() == 1 and
        grain_representations.at(2).getOrderParameterId() == 0)
      {
        result = true;
      }

    snprintf(buffer,
             sizeof(buffer),
             "Subtest 2 result for 'SimplifiedGrainManipulator::reassignGrains': %u\n",
             result);
    std::cout << buffer;

    pass = pass & result;
  }

  snprintf(buffer,
           sizeof(buffer),
           "Test result for 'SimplifiedGrainManipulator::reassignGrains': %u\n",
           pass);
  std::cout << buffer;

  return pass;
}
