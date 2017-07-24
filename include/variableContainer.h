// This class permits the access of a subset of indexed fields and gives an error if any non-allowed fields are requested
#ifndef VARIBLECONTAINER_H
#define VARIBLECONTAINER_H

#include "userInputParameters.h"

template <int dim, int degree, typename T>
class variableContainer
{
public:
    #include "typeDefs.h"

    // Constructors

    // Standard contructor, used for most situations
    variableContainer(const dealii::MatrixFree<dim,double> &data, std::vector<variable_info> _varInfoList);
    // Nonstandard constructor, used when only one index of "data" should be used, use with care!
    variableContainer(const dealii::MatrixFree<dim,double> &data, std::vector<variable_info> _varInfoList, unsigned int fixed_index);

    // Methods to get the value/grad/hess in the residual method (this is how the user gets these values in equations.h)
    T get_scalar_value(unsigned int global_variable_index) const;
    dealii::Tensor<1, dim, T > get_scalar_gradient(unsigned int global_variable_index) const;
    dealii::Tensor<2, dim, T > get_scalar_hessian(unsigned int global_variable_index) const;
    dealii::Tensor<1, dim, T > get_vector_value(unsigned int global_variable_index) const;
    dealii::Tensor<2, dim, T > get_vector_gradient(unsigned int global_variable_index) const;
    dealii::Tensor<3, dim, T > get_vector_hessian(unsigned int global_variable_index) const;

    // Methods to set the value residual and the gradient residual (this is how the user sets these values in equations.h)
    void set_scalar_value_residual_term(unsigned int global_variable_index, T val);
    void set_scalar_gradient_residual_term(unsigned int global_variable_index, dealii::Tensor<1, dim, T > grad);
    void set_vector_value_residual_term(unsigned int global_variable_index, dealii::Tensor<1, dim, T > val);
    void set_vector_gradient_residual_term(unsigned int global_variable_index, dealii::Tensor<2, dim, T > grad);

    // Initialize, read DOFs, and set evaulation flags for each variable
    void reinit_and_eval(const std::vector<vectorType*> &src, unsigned int cell);
    void reinit_and_eval_LHS(const vectorType &src, const std::vector<vectorType*> solutionSet, unsigned int cell, unsigned int var_being_solved);

    // Only initialize the FEEvaluation object for each variable (used for post-processing)
    void reinit(unsigned int cell);

    // Integrate the residuals and distribute from local to global
    void integrate_and_distribute(std::vector<vectorType*> &dst);
    void integrate_and_distribute_LHS(vectorType &dst, unsigned int var_being_solved);

    // The quadrature point index, a method to get the number of quadrature points per cell, and a method to get the xyz coordinates for the quadrature point
    unsigned int q_point;
    unsigned int get_num_q_points();
    dealii::Point<dim,T> get_q_point_location();

    // Method to obtain JxW (the weighted Jacobian)
    void get_JxW(dealii::AlignedVector<T> & JxW);

private:
    // The number of variables
    unsigned int num_var;

    // Vectors of the actual FEEvaluation objects for each active variable, split into scalar variables and vector variables for type reasons
    std::vector<dealii::FEEvaluation<dim,degree,degree+1,1,double> > scalar_vars;
    std::vector<dealii::FEEvaluation<dim,degree,degree+1,dim,double> > vector_vars;

    // Object containing some information about each variable (indices, whether the val/grad/hess is needed, etc)
    std::vector<variable_info> varInfoList;
};

#endif
