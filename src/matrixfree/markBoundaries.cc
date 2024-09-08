// methods to mark boundaries

#include "../../include/initialConditions.h"
#include "../../include/matrixFreePDE.h"

// methods to mark boundaries
// methods to mark boundaries
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::markBoundaries(
  parallel::distributed::Triangulation<dim> &tria) const
{
  for (const auto &cell : tria.active_cell_iterators())
    {
      // Mark all of the faces
      for (unsigned int face_number = 0; face_number < GeometryInfo<dim>::faces_per_cell;
           ++face_number)
        {
          for (unsigned int i = 0; i < dim; i++)
            {
              if (std::fabs(cell->face(face_number)->center()(i) - (0)) < 1e-12)
                {
                  cell->face(face_number)->set_boundary_id(2 * i);
                }
              else if (std::fabs(cell->face(face_number)->center()(i) -
                                 (userInputs.domain_size[i])) < 1e-12)
                {
                  cell->face(face_number)->set_boundary_id(2 * i + 1);
                }
            }
        }
    }
}

#include "../../include/matrixFreePDE_template_instantiations.h"
