// This class permits the access of a subset of indexed fields and gives an error if any non-allowed fields are requested
#ifndef VARIBLECONTAINER_H
#define VARIBLECONTAINER_H

#include "userInputParameters.h"

template <int dim, int degree, typename T>
class variableContainer
{
public:
    #include "typeDefs.h"

    // Constructor
    variableContainer(const dealii::MatrixFree<dim,double> &data, std::vector<variable_info> _varInfoList);

    // Methods to get the value/grad/hess in the residual method (this is how the user gets these values in equations.h)
    T get_scalar_value(unsigned int global_variable_index) const;
    dealii::Tensor<1, dim, T > get_scalar_gradient(unsigned int global_variable_index) const;
    dealii::Tensor<2, dim, T > get_scalar_hessian(unsigned int global_variable_index) const;
    dealii::Tensor<1, dim, T > get_vector_value(unsigned int global_variable_index) const;
    dealii::Tensor<2, dim, T > get_vector_gradient(unsigned int global_variable_index) const;
    dealii::Tensor<3, dim, T > get_vector_hessian(unsigned int global_variable_index) const;

    // Methods to set the value residual and the gradient residual (this is how the user sets these values in equations.h)
    void set_scalar_value_residual(unsigned int global_variable_index, T val);
    void set_scalar_gradient_residual(unsigned int global_variable_index, dealii::Tensor<1, dim, T > grad);
    void set_vector_value_residual(unsigned int global_variable_index, dealii::Tensor<1, dim, T > val);
    void set_vector_gradient_residual(unsigned int global_variable_index, dealii::Tensor<2, dim, T > grad);

    // Initialize, read DOFs, and set evaulation flags for each variable
    void reinit_and_eval(const std::vector<vectorType*> &src, unsigned int cell);
    void reinit_and_eval_LHS(const vectorType &src, const std::vector<vectorType*> solutionSet, unsigned int cell, unsigned int var_being_solved);

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

template <int dim, int degree, typename T>
void variableContainer<dim,degree,T>::get_JxW(dealii::AlignedVector<T> & JxW){
    if (scalar_vars.size() > 0){
        scalar_vars[0].fill_JxW_values(JxW);
    }
    else {
        vector_vars[0].fill_JxW_values(JxW);
    }
}

template <int dim, int degree, typename T>
unsigned int variableContainer<dim,degree,T>::get_num_q_points(){
    if (scalar_vars.size() > 0){
        return scalar_vars[0].n_q_points;
    }
    else {
        return vector_vars[0].n_q_points;
    }
}

template <int dim, int degree, typename T>
dealii::Point<dim, T> variableContainer<dim,degree,T>::get_q_point_location(){
    if (scalar_vars.size() > 0){
        return scalar_vars[0].quadrature_point(q_point);
    }
    else {
        return vector_vars[0].quadrature_point(q_point);
    }
}

template <int dim, int degree, typename T>
variableContainer<dim,degree,T>::variableContainer(const dealii::MatrixFree<dim,double> &data, std::vector<variable_info> _varInfoList)
{
    varInfoList = _varInfoList;

    num_var = varInfoList.size();

    for (unsigned int i=0; i < num_var; i++){
        if (varInfoList[i].var_needed){
            if (varInfoList[i].is_scalar){
                dealii::FEEvaluation<dim,degree,degree+1,1,double> var(data, i);
                scalar_vars.push_back(var);
            }
            else {
                dealii::FEEvaluation<dim,degree,degree+1,dim,double> var(data, i);
                vector_vars.push_back(var);
            }
        }
    }

}

template <int dim, int degree, typename T>
void variableContainer<dim,degree,T>::reinit_and_eval(const std::vector<vectorType*> &src, unsigned int cell){

    for (unsigned int i=0; i<num_var; i++){
        if (varInfoList[i].var_needed){
            if (varInfoList[i].is_scalar) {
                scalar_vars[varInfoList[i].scalar_or_vector_index].reinit(cell);
                scalar_vars[varInfoList[i].scalar_or_vector_index].read_dof_values(*src[i]);
                scalar_vars[varInfoList[i].scalar_or_vector_index].evaluate(varInfoList[i].need_value, varInfoList[i].need_gradient, varInfoList[i].need_hessian);
            }
            else {
                vector_vars[varInfoList[i].scalar_or_vector_index].reinit(cell);
                vector_vars[varInfoList[i].scalar_or_vector_index].read_dof_values(*src[i]);
                vector_vars[varInfoList[i].scalar_or_vector_index].evaluate(varInfoList[i].need_value, varInfoList[i].need_gradient, varInfoList[i].need_hessian);
            }
        }
    }
}

template <int dim, int degree, typename T>
void variableContainer<dim,degree,T>::reinit_and_eval_LHS(const vectorType &src, const std::vector<vectorType*> solutionSet, unsigned int cell, unsigned int var_being_solved){

    for (unsigned int i=0; i<num_var; i++){
        if (varInfoList[i].var_needed){
            if (varInfoList[i].is_scalar) {
                scalar_vars[varInfoList[i].scalar_or_vector_index].reinit(cell);
                if (i == var_being_solved ){
                    scalar_vars[varInfoList[i].scalar_or_vector_index].read_dof_values(src);
                }
                else{
                    scalar_vars[varInfoList[i].scalar_or_vector_index].read_dof_values(*solutionSet[i]);
                }
                scalar_vars[varInfoList[i].scalar_or_vector_index].evaluate(varInfoList[i].need_value, varInfoList[i].need_gradient, varInfoList[i].need_hessian);
            }
            else {
                vector_vars[varInfoList[i].scalar_or_vector_index].reinit(cell);
                if (i == var_being_solved){
                    vector_vars[varInfoList[i].scalar_or_vector_index].read_dof_values(src);
                }
                else {
                    vector_vars[varInfoList[i].scalar_or_vector_index].read_dof_values(*solutionSet[i]);
                }
                vector_vars[varInfoList[i].scalar_or_vector_index].evaluate(varInfoList[i].need_value, varInfoList[i].need_gradient, varInfoList[i].need_hessian);
            }
        }
    }

}

template <int dim, int degree, typename T>
void variableContainer<dim,degree,T>::integrate_and_distribute(std::vector<vectorType*> &dst){

    for (unsigned int i=0; i<num_var; i++){
        if (varInfoList[i].is_scalar) {
            scalar_vars[varInfoList[i].scalar_or_vector_index].integrate(varInfoList[i].value_residual, varInfoList[i].gradient_residual);
            scalar_vars[varInfoList[i].scalar_or_vector_index].distribute_local_to_global(*dst[i]);
        }
        else {
            vector_vars[varInfoList[i].scalar_or_vector_index].integrate(varInfoList[i].value_residual, varInfoList[i].gradient_residual);
            vector_vars[varInfoList[i].scalar_or_vector_index].distribute_local_to_global(*dst[i]);
        }
    }

}

template <int dim, int degree, typename T>
void variableContainer<dim,degree,T>::integrate_and_distribute_LHS(vectorType &dst, unsigned int var_being_solved){

    //integrate
    if (varInfoList[var_being_solved].is_scalar) {
    	scalar_vars[varInfoList[var_being_solved].scalar_or_vector_index].integrate(varInfoList[var_being_solved].value_residual, varInfoList[var_being_solved].gradient_residual);
    	scalar_vars[varInfoList[var_being_solved].scalar_or_vector_index].distribute_local_to_global(dst);
    }
    else {
    	vector_vars[varInfoList[var_being_solved].scalar_or_vector_index].integrate(varInfoList[var_being_solved].value_residual, varInfoList[var_being_solved].gradient_residual);
    	vector_vars[varInfoList[var_being_solved].scalar_or_vector_index].distribute_local_to_global(dst);
    }
}

// Need to add index checking to these functions so that an error is thrown if the index wasn't set
template <int dim, int degree, typename T>
T variableContainer<dim,degree,T>::get_scalar_value(unsigned int global_variable_index) const
{
    if (varInfoList[global_variable_index].need_value){
        return scalar_vars[varInfoList[global_variable_index].scalar_or_vector_index].get_value(q_point);
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<1, dim, T > variableContainer<dim,degree,T>::get_scalar_gradient(unsigned int global_variable_index) const
{
    if (varInfoList[global_variable_index].need_gradient){
        return scalar_vars[varInfoList[global_variable_index].scalar_or_vector_index].get_gradient(q_point);
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<2, dim, T > variableContainer<dim,degree,T>::get_scalar_hessian(unsigned int global_variable_index) const
{
    if (varInfoList[global_variable_index].need_hessian){
        return scalar_vars[varInfoList[global_variable_index].scalar_or_vector_index].get_hessian(q_point);
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<1, dim, T > variableContainer<dim,degree,T>::get_vector_value(unsigned int global_variable_index) const
{
    if (varInfoList[global_variable_index].need_value){
        return vector_vars[varInfoList[global_variable_index].scalar_or_vector_index].get_value(q_point);
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<2, dim, T > variableContainer<dim,degree,T>::get_vector_gradient(unsigned int global_variable_index) const
{
    if (varInfoList[global_variable_index].need_gradient){
        return vector_vars[varInfoList[global_variable_index].scalar_or_vector_index].get_gradient(q_point);
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<3, dim, T > variableContainer<dim,degree,T>::get_vector_hessian(unsigned int global_variable_index) const
{
    if (varInfoList[global_variable_index].need_hessian){
        return vector_vars[varInfoList[global_variable_index].scalar_or_vector_index].get_hessian(q_point);
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

template <int dim, int degree, typename T>
void variableContainer<dim,degree,T>::set_scalar_value_residual(unsigned int global_variable_index, T val){
    scalar_vars[varInfoList[global_variable_index].scalar_or_vector_index].submit_value(val,q_point);
}

template <int dim, int degree, typename T>
void variableContainer<dim,degree,T>::set_scalar_gradient_residual(unsigned int global_variable_index, dealii::Tensor<1, dim, T > grad){
    scalar_vars[varInfoList[global_variable_index].scalar_or_vector_index].submit_gradient(grad,q_point);
}

template <int dim, int degree, typename T>
void variableContainer<dim,degree,T>::set_vector_value_residual(unsigned int global_variable_index, dealii::Tensor<1, dim, T > val){
    vector_vars[varInfoList[global_variable_index].scalar_or_vector_index].submit_value(val,q_point);
}
template <int dim, int degree, typename T>
void variableContainer<dim,degree,T>::set_vector_gradient_residual(unsigned int global_variable_index, dealii::Tensor<2, dim, T > grad){
    vector_vars[varInfoList[global_variable_index].scalar_or_vector_index].submit_gradient(grad,q_point);
}

#endif
