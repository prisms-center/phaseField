#ifndef triangulation_handler_h
#define triangulation_handler_h

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria.h>

#include <prismspf/config.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <memory>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class handlers the generation and manipulation of triangulations.
 */
template <int dim>
class triangulationHandler
{
public:
  using Triangulation =
    std::conditional_t<dim == 1,
                       dealii::Triangulation<dim>,
                       dealii::parallel::distributed::Triangulation<dim>>;

  /**
   * \brief Constructor.
   */
  triangulationHandler(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Getter function for triangulation (constant reference).
   */
  [[nodiscard]] const Triangulation &
  get_triangulation() const;

  /**
   * \brief Return the global maximum level of the triangulation.
   */
  [[nodiscard]] unsigned int
  get_n_global_levels() const;

  /**
   * \brief Generate mesh.
   */
  void
  generate_mesh();

  /**
   * \brief Export triangulation to vtk.
   */
  void
  export_triangulation_as_vtk(const std::string &filename) const;

private:
  /**
   * \brief Mark the domain ids on the triangulation to get the proper mapping of
   * specified boundary conditions.
   */
  void
  mark_boundaries() const;

  /**
   * \brief Mark certain faces of the triangulation periodic.
   */
  void
  mark_periodic();

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief Triangulation.
   */
  std::unique_ptr<Triangulation> triangulation;
};

PRISMS_PF_END_NAMESPACE

#endif