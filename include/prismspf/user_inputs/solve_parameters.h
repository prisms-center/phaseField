#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_selector.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/config.h>

#include <execution>
#include <map>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that stores relevant linear solve information of a certain solve block.
 */
struct LinearSolverParameters : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int n_subsections = Numbers::default_subsections) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              n_subsections = Numbers::default_subsections) override;

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override;

  // Solver type. richardson|cg|bicgstab|gmres|fgmres|minres
  std::string solver_type = "cg";

  // Solver tolerance
  double tolerance = 1.0e-10;

  // Solver tolerance type
  SolverToleranceType tolerance_type = SolverToleranceType::RMSEPerField;

  // Max number of iterations for the linear solve
  unsigned int max_iterations = 100;

  // The multigrid depth
  unsigned int mg_depth = 1;

  // Preconditioner
  PreconditionerType preconditioner = PreconditionerType::None;

  // Solver AdditionalData structures
  dealii::PreconditionChebyshev<>::AdditionalData chebyshev_parameters;

  void
  declare_chebyshev(dealii::ParameterHandler &parameter_handler) const;

  void
  assign_chebyshev(dealii::ParameterHandler &parameter_handler);

  dealii::SolverRichardson<>::AdditionalData richardson_parameters;

  void
  declare_richardson(dealii::ParameterHandler &parameter_handler) const;

  void
  assign_richardson(dealii::ParameterHandler &parameter_handler);

  dealii::SolverBicgstab<>::AdditionalData bicgstab_parameters;

  void
  declare_bicgstab(dealii::ParameterHandler &parameter_handler) const;

  void
  assign_bicgstab(dealii::ParameterHandler &parameter_handler);

  dealii::SolverGMRES<>::AdditionalData gmres_parameters;

  void
  declare_gmres(dealii::ParameterHandler &parameter_handler) const;

  void
  assign_gmres(dealii::ParameterHandler &parameter_handler);

  // MinRes and CG do not have additional parameters, so we do not need to store them
  // here.
  //   dealii::SolverMinRes<>::AdditionalData minres_parameters;
  //   dealii::SolverCG<>::AdditionalData     cg_parameters;
  // FGMRES additional parameters are a subset of GMRES parameters, so we can just use
  // gmres_parameters for both solvers
  //   dealii::SolverFGMRES<>::AdditionalData     fgmres_parameters;
};

/**
 * @brief Struct that holds linear solver parameters.
 */
struct LinearSolveParameters : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int n_subsections = Numbers::default_subsections) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              n_subsections = Numbers::default_subsections) override;

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override;

  // Map of linear solve parameters for fields that require them
  std::map<unsigned int, LinearSolverParameters> linear_solvers;
};

/**
 * @brief Struct that stores relevant nonlinear solve information of a certain field
 */
struct NonlinearSolverParameters : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int n_subsections = Numbers::default_subsections) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              n_subsections = Numbers::default_subsections) override;

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override;

  // Nonlinear step length
  double step_length = 1.0;

  // Max number of iterations for the nonlinear solve
  unsigned int max_iterations = 100;

  // Tolerance value for the nonlinear solve
  double tolerance_value = 1.0e-10;
};

/**
 * @brief Struct that holds nonlinear solver parameters.
 */
struct NonlinearSolveParameters : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int n_subsections = Numbers::default_subsections) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              n_subsections = Numbers::default_subsections) override;

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override;

  // Map of nonlinear solve parameters for fields that require them
  std::map<Types::Index, NonlinearSolverParameters> nonlinear_solvers;
};

PRISMS_PF_END_NAMESPACE
