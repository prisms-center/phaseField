#ifndef matrix_free_operator_h
#define matrix_free_operator_h

#include <deal.II/base/vectorization.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/operators.h>

#include <core/type_enums.h>
#include <core/user_inputs/user_input_parameters.h>
#include <core/variable_attributes.h>
#include <core/variable_container.h>

/**
 * \brief This is the abstract base class for the matrix-free implementation of some PDE.
 *
 * \tparam dim The number of dimensions in the problem.
 * \tparam degree The polynomial degree of the shape functions.
 * \tparam number Datatype to use for `LinearAlgebra::distributed::Vector<number>`. Either
 * double or float.
 */
template <int dim, int degree, typename number>
class matrixFreeOperator
  : public dealii::MatrixFreeOperators::
      Base<dim, dealii::LinearAlgebra::distributed::Vector<number>>
{
public:
  /**
   * \brief Constructor.
   */
  matrixFreeOperator(const AttributesList &_attributes_list);

  /**
   * \brief Release all memory and return to state like having called the default
   * constructor.
   */
  void
  clear() override;

  /**
   * \brief Compute the residual of this operator. This is the b in Ax=b.
   */
  void
  compute_residual(dealii::LinearAlgebra::distributed::Vector<double>       &dst,
                   const dealii::LinearAlgebra::distributed::Vector<double> &src) const;

  /**
   * \brief Compute the diagonal of this operator.
   */
  void
  compute_diagonal() override;

protected:
  /**
   * \brief User-implemented class for the RHS of explicit equations.
   */
  virtual void
  compute_explicit_RHS(
    variableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc) const = 0;

  /**
   * \brief User-implemented class for the RHS of nonexplicit equations.
   */
  virtual void
  compute_nonexplicit_RHS(
    variableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc) const = 0;

  /**
   * \brief User-implemented class for the LHS of nonexplicit equations.
   */
  virtual void
  compute_nonexplicit_LHS(
    variableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc) const = 0;

private:
  /**
   * \brief Apply operator to src and add result in dst.
   */
  void
  apply_add(dealii::LinearAlgebra::distributed::Vector<number>       &dst,
            const dealii::LinearAlgebra::distributed::Vector<number> &src) const override;

  /**
   * \brief Local application of the operator.
   */
  void
  local_apply(const dealii::MatrixFree<dim, number>                    &data,
              dealii::LinearAlgebra::distributed::Vector<number>       &dst,
              const dealii::LinearAlgebra::distributed::Vector<number> &src,
              const std::pair<uint, uint>                              &cell_range) const;

  /**
   * \brief Local computation of the residual of the operator.
   */
  void
  compute_local_residual(const dealii::MatrixFree<dim, double>                    &data,
                         dealii::LinearAlgebra::distributed::Vector<double>       &dst,
                         const dealii::LinearAlgebra::distributed::Vector<double> &src,
                         const std::pair<uint, uint> &cell_range) const;

  /**
   * \brief Local computation of the diagonal of the operator.
   */
  void
  local_compute_diagonal(const dealii::MatrixFree<dim, number>              &data,
                         dealii::LinearAlgebra::distributed::Vector<number> &dst,
                         const uint                                         &dummy,
                         const std::pair<uint, uint> &cell_range) const;

  /**
   * \brief The attribute list of the relevant variables.
   */
  const AttributesList &attributes_list;
};

template <int dim, int degree, typename number>
matrixFreeOperator<dim, degree, number>::matrixFreeOperator(
  const AttributesList &_attributes_list)
  : dealii::MatrixFreeOperators::
      Base<dim, dealii::LinearAlgebra::distributed::Vector<number>>()
  , attributes_list(_attributes_list)
{}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::clear()
{
  dealii::MatrixFreeOperators::
    Base<dim, dealii::LinearAlgebra::distributed::Vector<number>>::clear();
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::local_apply(
  const dealii::MatrixFree<dim, number>                    &data,
  dealii::LinearAlgebra::distributed::Vector<number>       &dst,
  const dealii::LinearAlgebra::distributed::Vector<number> &src,
  const std::pair<uint, uint>                              &cell_range) const
{
  // Constructor for FEEvaluation objects
  variableContainer<dim, degree, number> variable_list(data,
                                                       attributes_list,
                                                       solveType::NONEXPLICIT_LHS);

  // Initialize, evaluate, and submit based on user function.
  variable_list.eval_local_operator(
    [this](variableContainer<dim, degree, number>                    &var_list,
           const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
    {
      this->compute_nonexplicit_LHS(var_list, q_point_loc);
    },
    dst,
    src,
    cell_range);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::apply_add(
  dealii::LinearAlgebra::distributed::Vector<number>       &dst,
  const dealii::LinearAlgebra::distributed::Vector<number> &src) const
{
  this->data->cell_loop(&matrixFreeOperator::local_apply, this, dst, src);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_residual(
  dealii::LinearAlgebra::distributed::Vector<double>       &dst,
  const dealii::LinearAlgebra::distributed::Vector<double> &src) const
{
  this->data->cell_loop(&matrixFreeOperator::compute_local_residual, this, dst, src);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_local_residual(
  const dealii::MatrixFree<dim, double>                    &data,
  dealii::LinearAlgebra::distributed::Vector<double>       &dst,
  const dealii::LinearAlgebra::distributed::Vector<double> &src,
  const std::pair<uint, uint>                              &cell_range) const
{
  // Constructor for FEEvaluation objects
  variableContainer<dim, degree, double> variable_list(data,
                                                       attributes_list,
                                                       solveType::NONEXPLICIT_RHS);

  // Initialize, evaluate, and submit based on user function.
  variable_list.eval_local_operator(
    [this](variableContainer<dim, degree, double>                    &var_list,
           const dealii::Point<dim, dealii::VectorizedArray<double>> &q_point_loc)
    {
      this->compute_nonexplicit_RHS(var_list, q_point_loc);
    },
    dst,
    src,
    cell_range);
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::compute_diagonal()
{
  this->inverse_diagonal_entries.reset(
    new DiagonalMatrix<dealii::LinearAlgebra::distributed::Vector<number>>());
  dealii::LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
    this->inverse_diagonal_entries->get_vector();
  this->data->initialize_dof_vector(inverse_diagonal);
  uint dummy = 0;
  this->data->cell_loop(&matrixFreeOperator::local_compute_diagonal,
                        this,
                        inverse_diagonal,
                        dummy);

  this->set_constrained_entries_to_one(inverse_diagonal);

  for (uint i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
    {
      Assert(inverse_diagonal.local_element(i) > 0.0,
             dealii::ExcMessage(
               "No diagonal entry in a positive definite operator should be zero"));
      inverse_diagonal.local_element(i) = 1.0 / inverse_diagonal.local_element(i);
    }
}

template <int dim, int degree, typename number>
void
matrixFreeOperator<dim, degree, number>::local_compute_diagonal(
  const dealii::MatrixFree<dim, number>              &data,
  dealii::LinearAlgebra::distributed::Vector<number> &dst,
  [[maybe_unused]] const uint                        &dummy,
  const std::pair<uint, uint>                        &cell_range) const
{
  // Field index
  uint field_index = 0;

  // Constructor for FEEvaluation objects
  variableContainer<dim, degree, number> variable_list(data,
                                                       attributes_list,
                                                       solveType::NONEXPLICIT_LHS);

  // Initialize, evaluate, and submit diagonal based on user function.
  variable_list.eval_local_diagonal(
    [this](variableContainer<dim, degree, number>                    &var_list,
           const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
    {
      this->compute_nonexplicit_LHS(var_list, q_point_loc);
    },
    dst,
    cell_range,
    field_index);
}

#endif