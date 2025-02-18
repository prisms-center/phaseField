#ifndef CUSTOM_PDE_H_
#define CUSTOM_PDE_H_

#include <prismspf/config.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/user_inputs/user_input_parameters.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This is a derived class of `matrixFreeOperator` where the user implements their
 * PDEs.
 *
 * \tparam dim The number of dimensions in the problem.
 * \tparam degree The polynomial degree of the shape functions.
 * \tparam number Datatype to use. Either double or float.
 */
template <int dim, int degree, typename number>
class customPDE : public matrixFreeOperator<dim, degree, number>
{
public:
  using scalarValue = dealii::VectorizedArray<number>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using scalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using vectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;

  /**
   * \brief Constructor.
   */
  customPDE(const userInputParameters<dim>                   &_user_inputs,
            const std::map<unsigned int, variableAttributes> &subset_attributes,
            const unsigned int                               &_current_index = 0)
    : matrixFreeOperator<dim, degree, number>(subset_attributes)
    , user_inputs(_user_inputs)
    , current_index(_current_index)
  {}

private:
  /**
   * \brief User-implemented class for the RHS of explicit equations.
   */
  void
  compute_explicit_RHS(variableContainer<dim, degree, number> &variable_list,
                       const dealii::Point<dim, dealii::VectorizedArray<number>>
                         &q_point_loc) const override;

  /**
   * \brief User-implemented class for the RHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_RHS(variableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, dealii::VectorizedArray<number>>
                            &q_point_loc) const override;

  /**
   * \brief User-implemented class for the LHS of nonexplicit equations.
   */
  void
  compute_nonexplicit_LHS(variableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, dealii::VectorizedArray<number>>
                            &q_point_loc) const override;

  /**
   * \brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  void
  compute_postprocess_explicit_RHS(
    variableContainer<dim, degree, number>                    &variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
    const override;

  /**
   * \brief The user-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief The current index that is being solved. For explicit equations, they are
   * solved concurrently, so this should not be used.
   */
  const unsigned int current_index;

  const static unsigned int          CIJ_tensor_size = (2 * dim) - 1 + (dim / 3);
  dealii::Tensor<2, CIJ_tensor_size> CIJ =
    user_inputs.user_constants.get_model_constant_elasticity_tensor("CIJ");
};

PRISMS_PF_END_NAMESPACE

#endif