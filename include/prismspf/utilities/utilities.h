// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>

#include <prismspf/config.h>

#include <array>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Positive moldulo (remainder)
 * returns the normal remainder. (c++ fmod is defined abnormally for negative numbers)
 */
template <typename Real, typename OtherReal>
inline Real
pmod(const Real &value, const OtherReal &modulus)
{
  using std::fmod;
  return fmod(fmod(value, modulus) + modulus, modulus);
}

/**
 * @brief Voigt notation index range
 */
template <unsigned int dim>
constexpr unsigned int voigt_tensor_size = (2 * dim) - 1 + (dim / 3);

/**
 * @brief Compute the stress with a given displacement and elasticity tensor. This assumes
 * that the provided parameters are in Voigt notation.
 */
template <unsigned int dim, typename T>
inline void
compute_stress(const dealii::Tensor<2, voigt_tensor_size<dim>, T> &elasticity_tensor,
               const dealii::Tensor<1, voigt_tensor_size<dim>, T> &strain,
               dealii::Tensor<1, voigt_tensor_size<dim>, T>       &stress)
{
  stress = elasticity_tensor * strain;
}

/**
 * @brief Compute the stress with a given displacement and elasticity tensor. Note: this
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
      const int xx_dir = 0;
      const int yy_dir = 1;
      const int zz_dir = 2;
      const int yz_dir = 3;
      const int xz_dir = 4;
      const int xy_dir = 5;

      epsilon[xx_dir] = strain[xx_dir][xx_dir];
      epsilon[yy_dir] = strain[yy_dir][yy_dir];
      epsilon[zz_dir] = strain[zz_dir][zz_dir];
      // In Voigt notation: epsilonngineering shear strain=zz_dir*strain
      epsilon[yz_dir] = strain[yy_dir][zz_dir] + strain[zz_dir][yy_dir];
      epsilon[xz_dir] = strain[xx_dir][zz_dir] + strain[zz_dir][xx_dir];
      epsilon[xy_dir] = strain[xx_dir][yy_dir] + strain[yy_dir][xx_dir];

      // Multiply elasticity_tensor and epsilon to get sigma
      sigma = elasticity_tensor * epsilon;

      stress[xx_dir][xx_dir] = sigma[xx_dir];
      stress[yy_dir][yy_dir] = sigma[yy_dir];
      stress[zz_dir][zz_dir] = sigma[zz_dir];

      stress[yy_dir][zz_dir] = sigma[yz_dir];
      stress[zz_dir][yy_dir] = sigma[yz_dir];

      stress[xx_dir][zz_dir] = sigma[xz_dir];
      stress[zz_dir][xx_dir] = sigma[xz_dir];

      stress[xx_dir][yy_dir] = sigma[xy_dir];
      stress[yy_dir][xx_dir] = sigma[xy_dir];
    }
  else if constexpr (dim == 2)
    {
      const int xx_dir = 0;
      const int yy_dir = 1;
      const int xy_dir = 2;

      epsilon[xx_dir] = strain[xx_dir][xx_dir];
      epsilon[yy_dir] = strain[yy_dir][yy_dir];
      // In Voigt notation: epsilonngineering shear strain=xy_dir*strain
      epsilon[xy_dir] = strain[xx_dir][yy_dir] + strain[yy_dir][xx_dir];

      // Multiply elasticity_tensor and epsilon to get sigma
      sigma = elasticity_tensor * epsilon;

      stress[xx_dir][xx_dir] = sigma[xx_dir];
      stress[yy_dir][yy_dir] = sigma[yy_dir];
      stress[xx_dir][yy_dir] = sigma[xy_dir];
      stress[yy_dir][xx_dir] = sigma[xy_dir];
    }
  else
    {
      const int xx_dir = 0;

      stress[xx_dir][xx_dir] = elasticity_tensor[xx_dir][xx_dir] * strain[xx_dir][xx_dir];
    }
}

/**
 * @brief Remove whitepace from strings
 */
inline std::string
strip_whitespace(const std::string &_text)
{
  std::string text = _text;
  boost::range::remove_erase_if(text, ::isspace);
  return text;
}

/**
 * @brief Convert bool to string.
 */
inline const char *
bool_to_string(bool boolean)
{
  return boolean ? "true" : "false";
}

/**
 * @brief Convert evaluation flags to string.
 */
inline std::string
eval_flags_to_string(dealii::EvaluationFlags::EvaluationFlags flag)
{
  std::string result;

  if ((flag & dealii::EvaluationFlags::EvaluationFlags::values) != 0U)
    {
      result += "values";
    }
  if ((flag & dealii::EvaluationFlags::EvaluationFlags::gradients) != 0U)
    {
      if (!result.empty())
        {
          result += " | ";
        }
      result += "gradients";
    }
  if ((flag & dealii::EvaluationFlags::EvaluationFlags::hessians) != 0U)
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

template <unsigned int dim, typename number>
inline std::vector<number>
dealii_point_to_vector(const dealii::Point<dim, number> &point)
{
  static_assert(dim < 4, "We only allow 3 space dimensions");

  std::vector<number> vec(3, 0.0);

  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-constant-array-index)

  for (unsigned int i = 0; i < dim; ++i)
    {
      vec[i] = point[i];
    }

  // NOLINTEND(cppcoreguidelines-pro-bounds-constant-array-index)

  return vec;
}

PRISMS_PF_END_NAMESPACE

#include <prismspf/utilities/vectorized_operations.h>