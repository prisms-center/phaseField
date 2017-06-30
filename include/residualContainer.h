// This class permits the access of a subset of indexed fields and gives an error if any non-allowed fields are requested
#ifndef RESIDUALCONTAINER_H
#define RESIDUALCONTAINER_H

#include "userInputParameters.h"

template <int dim, typename T>
class residualContainer
{
public:
    residualContainer(std::vector<bool> value_residual, std::vector<bool> gradient_residual, std::vector<fieldType> var_type);

    void set_scalar_value(unsigned int global_variable_index, T val);
    void set_scalar_gradient(unsigned int global_variable_index, dealii::Tensor<1, dim, T > grad);
    void set_vector_value(unsigned int global_variable_index, dealii::Tensor<1, dim, T > val);
    void set_vector_gradient(unsigned int global_variable_index, dealii::Tensor<2, dim, T > grad);

    T get_scalar_value(unsigned int global_variable_index) const;
    dealii::Tensor<1, dim, T > get_scalar_gradient(unsigned int global_variable_index) const;
    dealii::Tensor<1, dim, T > get_vector_value(unsigned int global_variable_index) const;
    dealii::Tensor<2, dim, T > get_vector_gradient(unsigned int global_variable_index) const;

private:
    // These vectors of indices convert from the global index to the index of each of the value/gradient/hessian vectors.
    std::vector<int> scalar_value_index;
    std::vector<int> scalar_gradient_index;
    std::vector<int> vector_value_index;
    std::vector<int> vector_gradient_index;

    std::vector<T> scalar_value;
    std::vector<dealii::Tensor<1, dim, T > > scalar_gradient;
    std::vector<dealii::Tensor<1, dim, T > > vector_value;
    std::vector<dealii::Tensor<2, dim, T > > vector_gradient;

};

template <int dim, typename T>
residualContainer<dim,T>::residualContainer(std::vector<bool> value_residual, std::vector<bool> gradient_residual, std::vector<fieldType> var_type)
{
    unsigned int num_scalar_values = 0;
    unsigned int num_scalar_gradients = 0;
    unsigned int num_vector_values = 0;
    unsigned int num_vector_gradients = 0;

    for (unsigned int var=0; var < value_residual.size(); var++){
        if (var_type[var] == SCALAR){
            if (value_residual[var]){
                scalar_value_index.push_back(num_scalar_values);
                num_scalar_values++;
            }
            else {
                scalar_value_index.push_back(-1);
            }
            if (gradient_residual[var]){
                scalar_gradient_index.push_back(num_scalar_gradients);
                num_scalar_gradients++;
            }
            else {
                scalar_gradient_index.push_back(-1);
            }
            vector_value_index.push_back(-1);
            vector_gradient_index.push_back(-1);
        }
        else {
            if (value_residual[var]){
                vector_value_index.push_back(num_vector_values);
                num_vector_values++;
            }
            else {
                vector_value_index.push_back(-1);
            }
            if (gradient_residual[var]){
                vector_gradient_index.push_back(num_vector_gradients);
                num_vector_gradients++;
            }
            else {
                vector_gradient_index.push_back(-1);
            }
            scalar_value_index.push_back(-1);
            scalar_gradient_index.push_back(-1);
        }
    }

    scalar_value.reserve(num_scalar_values);
    scalar_gradient.reserve(num_scalar_gradients);
    vector_value.reserve(num_vector_values);
    vector_gradient.reserve(num_vector_gradients);

}

template <int dim, typename T>
void residualContainer<dim,T>::set_scalar_value(unsigned int global_variable_index, T val){
    if (scalar_value_index[global_variable_index] != -1){
        scalar_value[scalar_value_index[global_variable_index]] = val;
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempt to set a residual value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a residual value is requested." << std::endl;
        abort();
    }
}

template <int dim, typename T>
void residualContainer<dim,T>::set_scalar_gradient(unsigned int global_variable_index, dealii::Tensor<1, dim, T > grad){
    if (scalar_gradient_index[global_variable_index] != -1){
        scalar_gradient[scalar_gradient_index[global_variable_index]] = grad;
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempt to set a residual value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a residual value is requested." << std::endl;
        abort();
    }
}
template <int dim, typename T>
void residualContainer<dim,T>::set_vector_value(unsigned int global_variable_index, dealii::Tensor<1, dim, T > val){
    if (vector_value_index[global_variable_index] != -1){
        vector_value[vector_value_index[global_variable_index]] = val;
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempt to set a residual value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a residual value is requested." << std::endl;
        abort();
    }
}
template <int dim, typename T>
void residualContainer<dim,T>::set_vector_gradient(unsigned int global_variable_index, dealii::Tensor<2, dim, T > grad){
    if (vector_gradient_index[global_variable_index] != -1){
        vector_gradient[vector_gradient_index[global_variable_index]] = grad;
    }
    else {
        std::cerr << "PRISMS-PF Error: Attempt to set a residual value that was not marked as needed in 'parameters.in'. Double-check the indices in user functions where a residual value is requested." << std::endl;
        abort();
    }
}


template <int dim, typename T>
T residualContainer<dim,T>::get_scalar_value(unsigned int global_variable_index) const
{
    return scalar_value[scalar_value_index[global_variable_index]];
}

template <int dim, typename T>
dealii::Tensor<1, dim, T > residualContainer<dim,T>::get_scalar_gradient(unsigned int global_variable_index) const
{
    return scalar_gradient[scalar_gradient_index[global_variable_index]];
}

template <int dim, typename T>
dealii::Tensor<1, dim, T > residualContainer<dim,T>::get_vector_value(unsigned int global_variable_index) const
{
    return vector_value[vector_value_index[global_variable_index]];
}

template <int dim, typename T>
dealii::Tensor<2, dim, T > residualContainer<dim,T>::get_vector_gradient(unsigned int global_variable_index) const
{
    return vector_gradient[vector_gradient_index[global_variable_index]];
}

#endif
