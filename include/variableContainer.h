// This class permits the access of a subset of indexed fields and gives an error if any non-allowed fields are requested
#ifndef VARIBLECONTAINER_H
#define VARIBLECONTAINER_H

#include "userInputParameters.h"

template <int dim, int degree, typename T>
class variableContainer
{
public:
    #include "typeDefs.h"

    variableContainer(const dealii::MatrixFree<dim,double> &data, std::vector<variable_info> _varInfoList);

    void set_scalar_value_residual(unsigned int global_variable_index, T val);
    void set_scalar_gradient_residual(unsigned int global_variable_index, dealii::Tensor<1, dim, T > grad);
    void set_vector_value_residual(unsigned int global_variable_index, dealii::Tensor<1, dim, T > val);
    void set_vector_gradient_residual(unsigned int global_variable_index, dealii::Tensor<2, dim, T > grad);

    T get_scalar_value(unsigned int global_variable_index) const;
    dealii::Tensor<1, dim, T > get_scalar_gradient(unsigned int global_variable_index) const;
    dealii::Tensor<2, dim, T > get_scalar_hessian(unsigned int global_variable_index) const;
    dealii::Tensor<1, dim, T > get_vector_value(unsigned int global_variable_index) const;
    dealii::Tensor<2, dim, T > get_vector_gradient(unsigned int global_variable_index) const;
    dealii::Tensor<3, dim, T > get_vector_hessian(unsigned int global_variable_index) const;

    unsigned int q_point;

    void reinit_and_eval(const std::vector<vectorType*> &src, unsigned int cell);

    void integrate_and_distribute(std::vector<vectorType*> &dst);

    unsigned int get_num_q_points();

    unsigned int num_var;

    dealii::Point<dim,T> get_q_point_location();

private:

    std::vector<dealii::FEEvaluation<dim,degree,degree+1,1,double> > scalar_vars;
    std::vector<dealii::FEEvaluation<dim,degree,degree+1,dim,double> > vector_vars;

    std::vector<variable_info> varInfoList;
};

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

template <int dim, int degree, typename T>
void variableContainer<dim,degree,T>::reinit_and_eval(const std::vector<vectorType*> &src, unsigned int cell){

    for (unsigned int i=0; i<num_var; i++){
        if (varInfoList[i].is_scalar) {
            scalar_vars[varInfoList[i].scalar_or_vector_index].reinit(cell);
            scalar_vars[varInfoList[i].scalar_or_vector_index].read_dof_values(*src[varInfoList[i].global_var_index]);
            scalar_vars[varInfoList[i].scalar_or_vector_index].evaluate(varInfoList[i].need_value, varInfoList[i].need_gradient, varInfoList[i].need_hessian);
        }
        else {
            vector_vars[varInfoList[i].scalar_or_vector_index].reinit(cell);
            vector_vars[varInfoList[i].scalar_or_vector_index].read_dof_values(*src[varInfoList[i].global_var_index]);
            vector_vars[varInfoList[i].scalar_or_vector_index].evaluate(varInfoList[i].need_value, varInfoList[i].need_gradient, varInfoList[i].need_hessian);
        }
    }

}

template <int dim, int degree, typename T>
void variableContainer<dim,degree,T>::integrate_and_distribute(std::vector<vectorType*> &dst){

    for (unsigned int i=0; i<num_var; i++){
        if (varInfoList[i].is_scalar) {
            scalar_vars[varInfoList[i].scalar_or_vector_index].integrate(varInfoList[i].value_residual, varInfoList[i].gradient_residual);
            scalar_vars[varInfoList[i].scalar_or_vector_index].distribute_local_to_global(*dst[varInfoList[i].global_var_index]);
        }
        else {
            vector_vars[varInfoList[i].scalar_or_vector_index].integrate(varInfoList[i].value_residual, varInfoList[i].gradient_residual);
            vector_vars[varInfoList[i].scalar_or_vector_index].distribute_local_to_global(*dst[varInfoList[i].global_var_index]);
        }
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
