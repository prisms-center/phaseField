// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/phase_field_tools.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/triangulation_handler.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class provides the context for the grid refiner with ptrs to all the
 * relevant dependencies.
 *
 * The context in this case refers to all the finite element machinery needed to refine
 * the grid and solutions. For example, it contains the triangulation handler and the DoF
 * handler.
 */
template <unsigned int dim, unsigned int degree, typename number>
class GridRefinementContext
{
public:
  /**
   * @brief Constructor.
   */
  GridRefinementContext(
    const UserInputParameters<dim>                         &_user_inputs,
    PhaseFieldTools<dim>                                   &_pf_tools,
    TriangulationHandler<dim>                              &_triangulation_handler,
    ConstraintHandler<dim, degree, number>                 &_constraint_handler,
    MatrixFreeContainer<dim, number>                       &_matrix_free_container,
    InvmHandler<dim, degree, number>                       &_invm_handler,
    SolutionHandler<dim, number>                           &_solution_handler,
    DofHandler<dim>                                        &_dof_handler,
    std::map<FieldInfo::TensorRank, dealii::FESystem<dim>> &_fe_system,
    const dealii::MappingQ1<dim>                           &_mapping,
    ElementVolumeContainer<dim, degree, number>            &_element_volume_container,
    const MGInfo<dim>                                      &_mg_info)
    : user_inputs(&_user_inputs)
    , pf_tools(&_pf_tools)
    , triangulation_handler(&_triangulation_handler)
    , constraint_handler(&_constraint_handler)
    , matrix_free_container(&_matrix_free_container)
    , invm_handler(&_invm_handler)
    , solution_handler(&_solution_handler)
    , dof_handler(&_dof_handler)
    , fe_system(&_fe_system)
    , mapping(&_mapping)
    , element_volume_container(&_element_volume_container)
    , mg_info(&_mg_info) {};

  /**
   * @brief Destructor.
   */
  ~GridRefinementContext() = default;

  /**
   * @brief Get the user-inputs.
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    Assert(user_inputs != nullptr, dealii::ExcNotInitialized());
    return *user_inputs;
  }

  [[nodiscard]] PhaseFieldTools<dim> &
  get_pf_tools() const
  {
    Assert(pf_tools != nullptr, dealii::ExcNotInitialized());
    return *pf_tools;
  }

  /**
   * @brief Get the constraint handler.
   */
  [[nodiscard]] ConstraintHandler<dim, degree, number> &
  get_constraint_handler() const
  {
    Assert(constraint_handler != nullptr, dealii::ExcNotInitialized());
    return *constraint_handler;
  }

  /**
   * @brief Get the matrix-free container.
   */
  [[nodiscard]] MatrixFreeContainer<dim, number> &
  get_matrix_free_container() const
  {
    Assert(matrix_free_container != nullptr, dealii::ExcNotInitialized());
    return *matrix_free_container;
  }

  /**
   * @brief Get the triangulation handler.
   */
  [[nodiscard]] TriangulationHandler<dim> &
  get_triangulation_handler() const
  {
    Assert(triangulation_handler != nullptr, dealii::ExcNotInitialized());
    return *triangulation_handler;
  }

  /**
   * @brief Get the invm handler.
   */
  [[nodiscard]] InvmHandler<dim, degree, number> &
  get_invm_handler() const
  {
    Assert(invm_handler != nullptr, dealii::ExcNotInitialized());
    return *invm_handler;
  }

  /**
   * @brief Get the dof handler.
   */
  [[nodiscard]] DofHandler<dim> &
  get_dof_handler() const
  {
    Assert(dof_handler != nullptr, dealii::ExcNotInitialized());
    return *dof_handler;
  }

  /**
   * @brief Get the mapping.
   */
  [[nodiscard]] const dealii::MappingQ1<dim> &
  get_mapping() const
  {
    Assert(mapping != nullptr, dealii::ExcNotInitialized());
    return *mapping;
  }

  /**
   * @brief Get the solution handler.
   */
  [[nodiscard]] SolutionHandler<dim, number> &
  get_solution_handler() const
  {
    Assert(solution_handler != nullptr, dealii::ExcNotInitialized());
    return *solution_handler;
  }

  /**
   * @brief Get the collection of finite element systems.
   */
  [[nodiscard]] const std::map<FieldInfo::TensorRank, dealii::FESystem<dim>> &
  get_finite_element_systems() const
  {
    Assert(fe_system != nullptr, dealii::ExcNotInitialized());
    Assert(!fe_system->empty(),
           dealii::ExcMessage("The finite element system map is empty."));
    return *fe_system;
  }

  /**
   * @brief Get the element volume container.
   */
  [[nodiscard]] ElementVolumeContainer<dim, degree, number> &
  get_element_volume_container() const
  {
    Assert(element_volume_container != nullptr, dealii::ExcNotInitialized());
    return *element_volume_container;
  }

  /**
   * @brief Get the multigrid info.
   */
  [[nodiscard]] const MGInfo<dim> &
  get_multigrid_info() const
  {
    Assert(mg_info != nullptr, dealii::ExcNotInitialized());
    return *mg_info;
  }

private:
  /**
   * @brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * @brief Phase field tools.
   */
  PhaseFieldTools<dim> *pf_tools;

  /**
   * @brief Triangulation handler.
   */
  TriangulationHandler<dim> *triangulation_handler;

  /**
   * @brief Constraint handler.
   */
  ConstraintHandler<dim, degree, number> *constraint_handler;

  /**
   * @brief Matrix-free object container.
   */
  MatrixFreeContainer<dim, number> *matrix_free_container;

  /**
   * @brief invm handler.
   */
  InvmHandler<dim, degree, number> *invm_handler;

  /**
   * @brief Solution handler.
   */
  SolutionHandler<dim, number> *solution_handler;

  /**
   * @brief DoF handler.
   */
  DofHandler<dim> *dof_handler;

  /**
   * @brief Collection of finite element systems. This is just a collection of two
   * FESystem's: one for scalar fields and one for vector fields. For now they both use
   * FE_Q finite elements.
   */
  std::map<FieldInfo::TensorRank, dealii::FESystem<dim>> *fe_system;

  /**
   * @brief Mappings to and from reference cell.
   */
  const dealii::MappingQ1<dim> *mapping;

  /**
   * @brief Element volume container.
   */
  ElementVolumeContainer<dim, degree, number> *element_volume_container;

  /**
   * @brief Multigrid information
   */
  const MGInfo<dim> *mg_info;
};

PRISMS_PF_END_NAMESPACE