// This class permits the access of a subset of indexed fields and gives an error if any non-allowed fields are requested
#ifndef VARIBLECONTAINER_H
#define VARIBLECONTAINER_H

#include "userInputParameters.h"

template <int dim, typename T>
class variableContainer
{
public:
    variableContainer(std::vector<bool> need_value, std::vector<bool> need_gradient, std::vector<bool> need_hessian, std::vector<fieldType> var_type);

    void set_scalar_value(unsigned int global_variable_index, T val);
    void set_scalar_gradient(unsigned int global_variable_index, dealii::Tensor<1, dim, T > grad);
    void set_scalar_hessian(unsigned int global_variable_index, dealii::Tensor<2, dim, T > hess);
    void set_vector_value(unsigned int global_variable_index, dealii::Tensor<1, dim, T > val);
    void set_vector_gradient(unsigned int global_variable_index, dealii::Tensor<2, dim, T > grad);
    void set_vector_hessian(unsigned int global_variable_index, dealii::Tensor<3, dim, T > hess);

    T get_scalar_value(unsigned int global_variable_index) const;
    dealii::Tensor<1, dim, T > get_scalar_gradient(unsigned int global_variable_index) const;
    dealii::Tensor<2, dim, T > get_scalar_hessian(unsigned int global_variable_index) const;
    dealii::Tensor<1, dim, T > get_vector_value(unsigned int global_variable_index) const;
    dealii::Tensor<2, dim, T > get_vector_gradient(unsigned int global_variable_index) const;
    dealii::Tensor<3, dim, T > get_vector_hessian(unsigned int global_variable_index) const;

private:
    // These vectors of indices convert from the global index to the index of each of the value/gradient/hessian vectors.
    std::vector<int> scalar_value_index;
    std::vector<int> scalar_gradient_index;
    std::vector<int> scalar_hessian_index;
    std::vector<int> vector_value_index;
    std::vector<int> vector_gradient_index;
    std::vector<int> vector_hessian_index;

    std::vector<T> scalar_value;
    std::vector<dealii::Tensor<1, dim, T > > scalar_gradient;
    std::vector<dealii::Tensor<2, dim, T > > scalar_hessian;

    std::vector<dealii::Tensor<1, dim, T > > vector_value;
    std::vector<dealii::Tensor<2, dim, T > > vector_gradient;
    std::vector<dealii::Tensor<3, dim, T > > vector_hessian;

};

template <int dim, typename T>
variableContainer<dim,T>::variableContainer(std::vector<bool> need_value, std::vector<bool> need_gradient, std::vector<bool> need_hessian, std::vector<fieldType> var_type)
{
    int num_scalar_values = 0;
    int num_scalar_gradients = 0;
    int num_scalar_hessians = 0;
    int num_vector_values = 0;
    int num_vector_gradients = 0;
    int num_vector_hessians = 0;

    for (unsigned int var=0; var < need_value.size(); var++){
        if (var_type[var] == SCALAR){
            if (need_value[var]){
                scalar_value_index.push_back(num_scalar_values);
                num_scalar_values++;
            }
            else {
                scalar_value_index.push_back(-1);
            }
            if (need_gradient[var]){
                scalar_gradient_index.push_back(num_scalar_gradients);
                num_scalar_gradients++;
            }
            else {
                scalar_gradient_index.push_back(-1);
            }
            if (need_hessian[var]){
                scalar_hessian_index.push_back(num_scalar_hessians);
                num_scalar_hessians++;
            }
            else {
                scalar_hessian_index.push_back(-1);
            }
            vector_value_index.push_back(-1);
            vector_gradient_index.push_back(-1);
            vector_hessian_index.push_back(-1);
        }
        else {
            if (need_value[var]){
                vector_value_index.push_back(num_vector_values);
                num_vector_values++;
            }
            else {
                vector_value_index.push_back(-1);
            }
            if (need_gradient[var]){
                vector_gradient_index.push_back(num_vector_gradients);
                num_vector_gradients++;
            }
            else {
                vector_gradient_index.push_back(-1);
            }
            if (need_hessian[var]){
                vector_hessian_index.push_back(num_vector_hessians);
                num_vector_hessians++;
            }
            else {
                vector_hessian_index.push_back(-1);
            }
            scalar_value_index.push_back(-1);
            scalar_gradient_index.push_back(-1);
            scalar_hessian_index.push_back(-1);
        }
    }

    scalar_value.reserve(num_scalar_values);
    scalar_gradient.reserve(num_scalar_gradients);
    scalar_hessian.reserve(num_scalar_hessians);
    vector_value.reserve(num_vector_values);
    vector_gradient.reserve(num_vector_gradients);
    vector_hessian.reserve(num_vector_hessians);

}

template <int dim, typename T>
void variableContainer<dim,T>::set_scalar_value(unsigned int global_variable_index, T val){
    scalar_value[scalar_value_index[global_variable_index]] = val;
}

template <int dim, typename T>
void variableContainer<dim,T>::set_scalar_gradient(unsigned int global_variable_index, dealii::Tensor<1, dim, T > grad){
    scalar_gradient[scalar_gradient_index[global_variable_index]] = grad;
}
template <int dim, typename T>
void variableContainer<dim,T>::set_scalar_hessian(unsigned int global_variable_index, dealii::Tensor<2, dim, T > hess){
    scalar_hessian[scalar_hessian_index[global_variable_index]] = hess;
}
template <int dim, typename T>
void variableContainer<dim,T>::set_vector_value(unsigned int global_variable_index, dealii::Tensor<1, dim, T > val){
    vector_value[vector_value_index[global_variable_index]] = val;
}
template <int dim, typename T>
void variableContainer<dim,T>::set_vector_gradient(unsigned int global_variable_index, dealii::Tensor<2, dim, T > grad){
    vector_gradient[vector_gradient_index[global_variable_index]] = grad;
}
template <int dim, typename T>
void variableContainer<dim,T>::set_vector_hessian(unsigned int global_variable_index, dealii::Tensor<3, dim, T > hess){
    vector_hessian[vector_hessian_index[global_variable_index]] = hess;
}

// Need to add index checking to these functions so that an error is thrown if the index wasn't set
template <int dim, typename T>
T variableContainer<dim,T>::get_scalar_value(unsigned int global_variable_index) const
{
    if (scalar_value_index[global_variable_index] != -1){
        return scalar_value[scalar_value_index[global_variable_index]];
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

template <int dim, typename T>
dealii::Tensor<1, dim, T > variableContainer<dim,T>::get_scalar_gradient(unsigned int global_variable_index) const
{
    if (scalar_gradient_index[global_variable_index] != -1){
        return scalar_gradient[scalar_gradient_index[global_variable_index]];
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

template <int dim, typename T>
dealii::Tensor<2, dim, T > variableContainer<dim,T>::get_scalar_hessian(unsigned int global_variable_index) const
{
    if (scalar_hessian_index[global_variable_index] != -1){
        return scalar_hessian[scalar_hessian_index[global_variable_index]];
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

template <int dim, typename T>
dealii::Tensor<1, dim, T > variableContainer<dim,T>::get_vector_value(unsigned int global_variable_index) const
{
    if (vector_value_index[global_variable_index] != -1){
        return vector_value[vector_value_index[global_variable_index]];
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

template <int dim, typename T>
dealii::Tensor<2, dim, T > variableContainer<dim,T>::get_vector_gradient(unsigned int global_variable_index) const
{
    if (vector_gradient_index[global_variable_index] != -1){
        return vector_gradient[vector_gradient_index[global_variable_index]];
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

template <int dim, typename T>
dealii::Tensor<3, dim, T > variableContainer<dim,T>::get_vector_hessian(unsigned int global_variable_index) const
{
    if (vector_hessian_index[global_variable_index] != -1){
        return vector_hessian[vector_hessian_index[global_variable_index]];
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a variable value is requested." << std::endl;
        abort();
    }
}

#endif
