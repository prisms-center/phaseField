// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef utilities_h
#define utilities_h

#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/config.h>

#include <array>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Convert given scalar to vectorized array.
 */
template <typename number>
constexpr auto
constV(number value)
{
  return dealii::make_vectorized_array<number>(value);
}

/**
 * \brief Convert given vector to a tensor of vectorized arrays.
 */
template <int dim, typename number>
constexpr auto
constT(const std::array<number, dim> &vector)
{
  dealii::Tensor<1, dim, dealii::VectorizedArray<number>> tensor;

  // Populate the Tensor with vectorized arrays
  for (std::size_t i = 0; i < dim; ++i)
    {
      tensor[i] = dealii::make_vectorized_array<number>(vector[i]);
    }

  return tensor;
}

/**
 * \brief Convert bool to string.
 */
inline const char *
bool_to_string(bool b)
{
  return b ? "true" : "false";
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

#endif