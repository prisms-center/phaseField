// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/affine_constraints.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/grid_refiner_criterion.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Internal enum for various triangulation types.
 */
enum TriangulationType : std::uint8_t
{
  Rectangular,
  Spherical,
  Custom
};

/**
 * @brief Triangulation base class.
 */
template <unsigned int dim>
struct Mesh
{
  using Triangulation =
    std::conditional_t<dim == 1,
                       dealii::Triangulation<dim>,
                       dealii::parallel::distributed::Triangulation<dim>>;

  /**
   * @brief Constructor.
   */
  Mesh() = default;

  /**
   * @brief Generate the mesh.
   */
  virtual void
  generate_mesh(Triangulation &triangulation) const = 0;

  /**
   * @brief Mark the boundaries of the mesh.
   */
  virtual void
  mark_boundaries(Triangulation &triangulation) const;

  /**
   * @brief Mark the periodic faces of the mesh.
   */
  void
  mark_periodic(Triangulation &triangulation);

  /**
   * @brief Mark the periodic faces of the mesh.
   */
  template <typename number>
  void
  mark_periodic(dealii::DoFHandler<dim>           &dof_handler,
                dealii::AffineConstraints<number> &constraints);

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  virtual void
  declare_parameters(dealii::ParameterHandler &parameter_handler) const = 0;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  virtual void
  assign_parameters(dealii::ParameterHandler &parameter_handler) = 0;

  /**
   * @brief Validate.
   */
  virtual void
  validate() const = 0;

  /**
   * @brief Periodicity set.
   *
   * This set contains all the information about periodicity for the system. Each entry
   * contain tuple with three entries: the first boundary id, the second boundary id,
   * and the direction.
   *
   * For example, (0,1,0) would link the 0th and 1st boundaries with periodicity in the
   * x-direction.
   */
  std::set<std::tuple<unsigned int, unsigned int, unsigned int>> periodicity_set;

  /**
   * @brief Triangulation periodicity vector.
   */
  std::vector<dealii::GridTools::PeriodicFacePair<typename Triangulation::cell_iterator>>
    triangulation_periodicity_vector;

  /**
   * @brief DoFHandler periodicity vector.
   */
  std::vector<
    dealii::GridTools::PeriodicFacePair<typename dealii::DoFHandler<dim>::cell_iterator>>
    dof_handler_periodicity_vector;
};

/**
 * @brief Class for rectangular mesh parameters.
 */
template <unsigned int dim>
struct RectangularMesh : public Mesh<dim>
{
  using Triangulation = typename Mesh<dim>::Triangulation;

  /**
   * @brief Constructor.
   */
  RectangularMesh() = default;

  /**
   * @brief Constructor.
   */
  RectangularMesh(dealii::Tensor<1, dim, double> _upper_bound,
                  dealii::Tensor<1, dim, double> _lower_bound,
                  std::array<unsigned int, dim>  _subdivisions);

  /**
   * @brief Generate the mesh.
   */
  void
  generate_mesh(Triangulation &triangulation) const override;

  /**
   * @brief Mark the boundaries of the mesh.
   */
  void
  mark_boundaries(Triangulation &triangulation) const override;

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler) const override;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler) override;

  /**
   * @brief Validate
   */
  void
  validate() const override;

  /**
   * @brief Upper bound point.
   */
  dealii::Tensor<1, dim, double> upper_bound;

  /**
   * @brief Lower bound point.
   */
  dealii::Tensor<1, dim, double> lower_bound;

  /**
   * @brief Mesh subdivisions in each cartesian direction.
   */
  std::vector<unsigned int> subdivisions = std::vector<unsigned int>(dim, 1);
};

/**
 * @brief Class for spherical mesh parameters.
 */
template <unsigned int dim>
struct SphericalMesh : public Mesh<dim>
{
  using Triangulation = typename Mesh<dim>::Triangulation;

  /**
   * @brief Constructor.
   */
  SphericalMesh() = default;

  /**
   * @brief Constructor.
   */
  explicit SphericalMesh(double _radius);

  /**
   * @brief Generate the mesh.
   */
  void
  generate_mesh(Triangulation &triangulation) const;

  /**
   * @brief Mark the boundaries of the mesh.
   */
  void
  mark_boundaries(Triangulation &triangulation) const override;

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler) const override;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler) override;

  /**
   * @brief Validate
   */
  void
  validate() const override;

  /**
   * @brief Radius of the spherical domain.
   */
  double radius = 0.0;
};

/**
 * @brief Struct that holds spatial discretization parameters.
 */
template <unsigned int dim>
struct SpatialDiscretization : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  predeclare(dealii::ParameterHandler &parameter_handler) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  preassign(dealii::ParameterHandler &parameter_handler) override;

  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) override;

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override;

  /**
   * @brief Whether the provided increment is a valid grid refinement step.
   */
  [[nodiscard]] bool
  should_refine_mesh(unsigned int increment) const;

  // Mesh object
  std::unique_ptr<Mesh<dim>> mesh;

  // Global refinement of mesh
  unsigned int global_refinement = 0;

  // Whether adaptive meshing (AMR) is enabled
  bool has_adaptivity = false;

  // Maximum global refinement for AMR
  unsigned int max_refinement = 0;

  // Minimum global refinement for AMR
  unsigned int min_refinement = 0;

  // The number of steps between remeshing
  unsigned int remeshing_period = UINT_MAX;

  // The criteria used for remeshing
  std::map<std::string, RefinementCriterion> refinement_criteria;
};

PRISMS_PF_END_NAMESPACE
