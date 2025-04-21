// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>

#include <prismspf/config.h>

#include <array>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Convert given scalar to vectorized array.
 */
template <typename number, typename T>
constexpr auto
constV(T value)
{
  return dealii::make_vectorized_array<number>(static_cast<number>(value));
}

/**
 * \brief Convert given vector to a tensor of vectorized arrays.
 */
template <unsigned int dim, typename number>
constexpr auto
constT(const std::array<number, dim> &vector)
{
  dealii::Tensor<1, dim, dealii::VectorizedArray<number>> tensor;

  // Populate the Tensor with vectorized arrays
  for (int i = 0; i < dim; ++i)
    {
      tensor[i] = dealii::make_vectorized_array<number>(vector[i]);
    }

  return tensor;
}

/**
 * \brief Voigt notation index range
 */
template <unsigned int dim>
constexpr unsigned int voigt_tensor_size = (2 * dim) - 1 + (dim / 3);

/**
 * \brief Compute the stress with a given displacement and elasticity tensor. This assumes
 * that the provided parameters are in Voigt notation.
 */
template <unsigned int dim, typename T>
inline void
compute_stress(const dealii::Tensor<2, voigt_tensor_size<dim>, T> &elasticity_tensor,
               const dealii::Tensor<1, voigt_tensor_size<dim>, T> &strain,
               dealii::Tensor<1, (2 * dim) - 1 + (dim / 3), T>    &stress)
{
  stress = elasticity_tensor * strain;
}

/**
 * \brief Compute the stress with a given displacement and elasticity tensor. Note: this
 * functions internally converts to Voigt notation.
 */
template <unsigned int dim, typename T>
inline void
compute_stress(const dealii::Tensor<2, voigt_tensor_size<dim>, T> &elasticity_tensor,
               const dealii::Tensor<2, dim, T>                    &strain,
               dealii::Tensor<2, dim, T>                          &stress)
{
  dealii::Tensor<1, voigt_tensor_size<dim>, T> sigma;
  dealii::Tensor<1, voigt_tensor_size<dim>, T> epsilon;

  if constexpr (dim == 3)
    {
      epsilon[0] = strain[0][0];
      epsilon[1] = strain[1][1];
      epsilon[2] = strain[2][2];
      // In Voigt notation: epsilonngineering shear strain=2*strain
      epsilon[3] = strain[1][2] + strain[2][1];
      epsilon[4] = strain[0][2] + strain[2][0];
      epsilon[5] = strain[0][1] + strain[1][0];

      // Multiply elasticity_tensor and epsilon to get sigma
      sigma = elasticity_tensor * epsilon;

      stress[0][0] = sigma[0];
      stress[1][1] = sigma[1];
      stress[2][2] = sigma[2];

      stress[1][2] = sigma[3];
      stress[2][1] = sigma[3];

      stress[0][2] = sigma[4];
      stress[2][0] = sigma[4];

      stress[0][1] = sigma[5];
      stress[1][0] = sigma[5];
    }
  else if constexpr (dim == 2)
    {
      epsilon[0] = strain[0][0];
      epsilon[1] = strain[1][1];
      // In Voigt notation: epsilonngineering shear strain=2*strain
      epsilon[2] = strain[0][1] + strain[1][0];

      // Multiply elasticity_tensor and epsilon to get sigma
      sigma = elasticity_tensor * epsilon;

      stress[0][0] = sigma[0];
      stress[1][1] = sigma[1];
      stress[0][1] = sigma[2];
      stress[1][0] = sigma[2];
    }
  else
    {
      stress[0][0] = elasticity_tensor[0][0] * strain[0][0];
    }
}

/**
 * \brief Compute the stress with a given displacement and elasticity tensor. Note: this
 * functions internally converts to Voigt notation.
 */
template <unsigned int dim, typename T>
inline void
compute_stress(
  const T (&elasticity_tensor)[voigt_tensor_size<dim>][voigt_tensor_size<dim>],
  const T (&strain)[][dim],
  T (&stress)[][dim])
{
  std::array<T, voigt_tensor_size<dim>> sigma {};
  std::array<T, voigt_tensor_size<dim>> epsilon {};

  if constexpr (dim == 3)
    {
      epsilon[0] = strain[0][0];
      epsilon[1] = strain[1][1];
      epsilon[2] = strain[2][2];
      // Voight notation to shear strain = 2 * strain
      epsilon[3] = strain[1][2] + strain[2][1];
      epsilon[4] = strain[0][2] + strain[2][0];
      epsilon[5] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < voigt_tensor_size<dim>; i++)
        {
          sigma[i] = 0.0;
          for (unsigned int j = 0; j < voigt_tensor_size<dim>; j++)
            {
              sigma[i] += elasticity_tensor[i][j] * epsilon[j];
            }
        }
      stress[0][0] = sigma[0];
      stress[1][1] = sigma[1];
      stress[2][2] = sigma[2];
      stress[1][2] = sigma[3];
      stress[0][2] = sigma[4];
      stress[0][1] = sigma[5];
      stress[2][1] = sigma[3];
      stress[2][0] = sigma[4];
      stress[1][0] = sigma[5];
    }
  else if constexpr (dim == 2)
    {
      epsilon[0] = strain[0][0];
      epsilon[1] = strain[1][1];
      // In Voigt notation: Engineering shear strain=2*strain
      epsilon[2] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < voigt_tensor_size<dim>; i++)
        {
          sigma[i] = 0.0;
          for (unsigned int j = 0; j < voigt_tensor_size<dim>; j++)
            {
              sigma[i] += elasticity_tensor[i][j] * epsilon[j];
            }
        }
      stress[0][0] = sigma[0];
      stress[1][1] = sigma[1];
      stress[0][1] = sigma[2];
      stress[1][0] = sigma[2];
    }
  else
    {
      epsilon[0]   = strain[0][0];
      sigma[0]     = elasticity_tensor[0][0] * epsilon[0];
      stress[0][0] = sigma[0];
    }
}

/**
 * \brief Compute the stress with a given displacement and elasticity tensor. Note: this
 * functions internally converts to Voigt notation.
 */
template <unsigned int dim, typename T>
inline void
compute_stress(const dealii::Tensor<2, voigt_tensor_size<dim>, T> &elasticity_tensor,
               const T (&strain)[][dim],
               T (&stress)[][dim])
{
  std::array<T, voigt_tensor_size<dim>> sigma {};
  std::array<T, voigt_tensor_size<dim>> epsilon {};

  if (dim == 3)
    {
      epsilon[0] = strain[0][0];
      epsilon[1] = strain[1][1];
      epsilon[2] = strain[2][2];
      // In Voigt notation: epsilonngineering shear strain=2*strain
      epsilon[3] = strain[1][2] + strain[2][1];
      epsilon[4] = strain[0][2] + strain[2][0];
      epsilon[5] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < voigt_tensor_size<dim>; i++)
        {
          sigma[i] = 0.0;
          for (unsigned int j = 0; j < voigt_tensor_size<dim>; j++)
            {
              sigma[i] += elasticity_tensor[i][j] * epsilon[j];
            }
        }
      stress[0][0] = sigma[0];
      stress[1][1] = sigma[1];
      stress[2][2] = sigma[2];
      stress[1][2] = sigma[3];
      stress[0][2] = sigma[4];
      stress[0][1] = sigma[5];
      stress[2][1] = sigma[3];
      stress[2][0] = sigma[4];
      stress[1][0] = sigma[5];
    }
  else if (dim == 2)
    {
      epsilon[0] = strain[0][0];
      epsilon[1] = strain[1][1];
      // In Voigt notation: epsilonngineering shear strain=2*strain
      epsilon[2] = strain[0][1] + strain[1][0];
      for (unsigned int i = 0; i < voigt_tensor_size<dim>; i++)
        {
          sigma[i] = 0.0;
          for (unsigned int j = 0; j < voigt_tensor_size<dim>; j++)
            {
              sigma[i] += elasticity_tensor[i][j] * epsilon[j];
            }
        }
      stress[0][0] = sigma[0];
      stress[1][1] = sigma[1];
      stress[0][1] = sigma[2];
      stress[1][0] = sigma[2];
    }
  else
    {
      epsilon[0]   = strain[0][0];
      sigma[0]     = elasticity_tensor[0][0] * epsilon[0];
      stress[0][0] = sigma[0];
    }
}

/**
 * \brief Remove whitepace from strings
 */
inline std::string
strip_whitespace(const std::string &_text)
{
  std::string text = _text;
  boost::range::remove_erase_if(text, ::isspace);
  return text;
}

/**
 * \brief Convert bool to string.
 */
inline const char *
bool_to_string(bool boolean)
{
  return boolean ? "true" : "false";
}

/**
 * \brief Convert evaluation flags to string.
 */
inline std::string
eval_flags_to_string(dealii::EvaluationFlags::EvaluationFlags flag)
{
  std::string result;

  if ((flag & dealii::EvaluationFlags::values) != 0U)
    {
      result += "values";
    }
  if ((flag & dealii::EvaluationFlags::gradients) != 0U)
    {
      if (!result.empty())
        {
          result += " | ";
        }
      result += "gradients";
    }
  if ((flag & dealii::EvaluationFlags::hessians) != 0U)
    {
      if (!result.empty())
        {
          result += " | ";
        }
      result += "hessians";
    }

  if (result.empty())
    {
      return "nothing";
    }

  return result;
}

PRISMS_PF_END_NAMESPACE
