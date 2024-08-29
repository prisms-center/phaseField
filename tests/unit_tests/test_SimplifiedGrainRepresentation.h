#include "../../include/SimplifiedGrainRepresentation.h"

template <int dim, typename T>
bool
unitTest<dim, T>::test_SimplifiedGrainRepresentation()
{
  char buffer[100];

  std::cout << "\nTesting 'SimplifiedGrainRepresentation'... " << std::endl;

  bool pass = true;

  // Subtest 1 (single element)
  {
    GrainSet<dim> test_grain_set;

    test_grain_set.setGrainIndex(2);

    std::vector<dealii::Point<dim>> vertex_set(dealii::Utilities::fixed_power<dim>(2.0));
    {
      dealii::Point<dim> p(0.0, 0.75);
      vertex_set[0] = p;
    }
    {
      dealii::Point<dim> p(0.25, 0.75);
      vertex_set[1] = p;
    }
    {
      dealii::Point<dim> p(0.0, 1.0);
      vertex_set[2] = p;
    }
    {
      dealii::Point<dim> p(0.25, 1.0);
      vertex_set[3] = p;
    }

    test_grain_set.addVertexList(vertex_set);

    SimplifiedGrainRepresentation<dim> simplified_grain_representation(test_grain_set);

    std::cout << "Centroid: " << simplified_grain_representation.getCenter() << std::endl;
    std::cout << "Radius: " << simplified_grain_representation.getRadius() << std::endl;

    bool result = false;
    if ((std::abs(simplified_grain_representation.getCenter()(0) - 0.125) < 1.0e-10) and
        (std::abs(simplified_grain_representation.getCenter()(1) - 0.875) < 1.0e-10) and
        (std::abs(simplified_grain_representation.getRadius() - 0.125 * std::sqrt(2.0)) <
         1.0e-10))
      {
        result = true;
      }

    snprintf(buffer,
             sizeof(buffer),
             "Subtest 1 result for 'SimplifiedGrainRepresentation': %u\n",
             result);
    std::cout << buffer;
    pass = pass and result;
  }

  // Subtest 2 (L-shaped grain)
  {
    GrainSet<dim> test_grain_set;

    test_grain_set.setGrainIndex(2);

    {
      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(0.0, 0.75);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(0.25, 0.75);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(0.0, 1.0);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(0.25, 1.0);
        vertex_set[3] = p;
      }

      test_grain_set.addVertexList(vertex_set);
    }
    {
      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(0.0, 0.5);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(0.25, 0.5);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(0.0, 0.75);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(0.25, 0.75);
        vertex_set[3] = p;
      }

      test_grain_set.addVertexList(vertex_set);
    }
    {
      std::vector<dealii::Point<dim>> vertex_set(
        dealii::Utilities::fixed_power<dim>(2.0));
      {
        dealii::Point<dim> p(0.25, 0.75);
        vertex_set[0] = p;
      }
      {
        dealii::Point<dim> p(0.5, 0.75);
        vertex_set[1] = p;
      }
      {
        dealii::Point<dim> p(0.25, 1.0);
        vertex_set[2] = p;
      }
      {
        dealii::Point<dim> p(0.5, 1.0);
        vertex_set[3] = p;
      }

      test_grain_set.addVertexList(vertex_set);
    }

    SimplifiedGrainRepresentation<dim> simplified_grain_representation(test_grain_set);

    std::cout << "Centroid: " << simplified_grain_representation.getCenter() << std::endl;
    std::cout << "Radius: " << simplified_grain_representation.getRadius() << std::endl;

    bool result = false;

    double centroid_x = (0.125 * 3.0 + 0.25) / 3.0;
    double centroid_y = (0.875 * 3.0 - 0.25) / 3.0;
    double radius     = std::sqrt(dealii::Utilities::fixed_power<2>(centroid_x - 0.5) +
                              dealii::Utilities::fixed_power<2>(centroid_y - 1.0));

    if ((std::abs(simplified_grain_representation.getCenter()(0) - centroid_x) <
         1.0e-10) and
        (std::abs(simplified_grain_representation.getCenter()(1) - centroid_y) <
         1.0e-10) and
        (std::abs(simplified_grain_representation.getRadius() - radius) < 1.0e-10))
      {
        result = true;
      }

    snprintf(buffer,
             sizeof(buffer),
             "Subtest 2 result for 'SimplifiedGrainRepresentation': %u\n",
             result);
    std::cout << buffer;
    pass = pass and result;
  }

  snprintf(buffer,
           sizeof(buffer),
           "Test result for 'SimplifiedGrainRepresentation': %u\n",
           pass);
  std::cout << buffer;

  return pass;
}
