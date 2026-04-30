// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>

#include <prismspf/utilities/vectorized_operations.h>

#include <prismspf/config.h>

#include <string>
#include <string_view>
#include <unordered_map>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali-manager.h>
#  include <caliper/cali.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class does the boilerplate that precedes and succeeds the simulation.
 *
 * @note In this class we initialize MPI, so everything that follows executes per process.
 */
class MPIInitFinalize : public dealii::Utilities::MPI::MPI_InitFinalize
{
public:
  MPIInitFinalize(int &argc, char **&argv, unsigned int _max_n_threads = 0);

  ~MPIInitFinalize();

private:
  /**
   * @brief Number of threads to try and use.
   *
   * @note If PRISMS_PF_THREADS=OFF this is always set to 1 unless otherwise specified.
   */
#ifdef PRISMS_PF_THREADS
  static constexpr unsigned int max_n_threads = dealii::numbers::invalid_unsigned_int;
#else
  static constexpr unsigned int max_n_threads = 1;
#endif

  /**
   * @brief Caliper config manager
   */
#ifdef PRISMS_PF_WITH_CALIPER
  cali::ConfigManager mgr;
#endif
};

/**
 * @brief Helper function that converts a string to some type, given a mapping.
 */
template <typename Type>
Type
string_to_type(const std::string                           &string,
               const std::unordered_map<std::string, Type> &table)
{
  auto iterator = table.find(string);
  AssertThrow(iterator != table.end(),
              dealii::ExcMessage("Unknown table entry: " + string));
  return iterator->second;
}

/**
 * @brief Helper function that converts a string to some type pair, given two mappings and
 * a set of delimiters.
 *
 * Note that this function requires the general form of one of the below cases.
 * Parentheses are used for delimiters for simplicity.
 *
 * Case 1:
 *  Input - Type(OtherType)
 *  Output - {OtherType, Type}
 * Case 2:
 *  Input - OtherType
 *  Output - {OtherType, Type}
 *
 * The input string would have the types replaced with their respective mappings.
 *
 */
template <typename Type, typename OtherType>
std::pair<OtherType, Type>
string_to_type_pair_with_delimiters(
  const std::string                                &string,
  const std::unordered_map<std::string, Type>      &table,
  const std::unordered_map<std::string, OtherType> &other_table,
  const std::pair<char, char>                      &delimiters = {'(', ')'})
{
  // Check to see if the string has any delimiters
  const auto opening_delimiter_position = string.find(delimiters.first);
  const auto closing_delimiter_position = string.find(delimiters.second);

  bool has_opening_delimiter = opening_delimiter_position != std::string_view::npos;
  bool has_closing_delimiter = closing_delimiter_position != std::string_view::npos;

  // Some checks for malformed input
  AssertThrow(closing_delimiter_position >= opening_delimiter_position,
              dealii::ExcMessage(
                "Opening delimiter must precede the closing delimiter. You had " +
                string));
  AssertThrow(has_closing_delimiter == has_opening_delimiter,
              dealii::ExcMessage("You must either have a delimiter pair or not. You seem "
                                 "to only have one delimiter. You had " +
                                 string));

  // Case 2 - no delimiters
  if (!has_closing_delimiter && !has_opening_delimiter)
    {
      const auto second_type = string_to_type("", table);
      const auto first_type  = string_to_type(string, other_table);
      return {first_type, second_type};
    }

  // More checks for malformed input
  AssertThrow(closing_delimiter_position == string.size() - 1,
              dealii::ExcMessage(
                "Closing delimiter must be at the end of the string. You had " + string));

  // Case 1 - delimiters
  std::string_view input_view {string};
  std::string_view outer_string = input_view.substr(0, opening_delimiter_position);
  std::string_view inner_string =
    input_view.substr(opening_delimiter_position + 1,
                      closing_delimiter_position - opening_delimiter_position - 1);

  AssertThrow(!outer_string.empty() && !inner_string.empty(),
              dealii::ExcMessage(
                "Inner and outer map entries must not be empty. You had " + string));
  // NOTE: This will allocate the string views to strings on the heap. Hopefully this
  // isn't so bad in performance, but I didn't wanna support heterogeneous lookup tables.
  // I can't imagine this becoming a performance problem right now since we only deal with
  // strings in the initialization of simulations; however, if it does, it shouldn't be
  // two hard to support heterogeneous lookup or a table that uses std::string_view.
  const auto second_type = string_to_type(std::string(outer_string), table);
  const auto first_type  = string_to_type(std::string(inner_string), other_table);
  return {first_type, second_type};
}

/**
 * @brief Positive moldulo (remainder)
 * returns the normal remainder. (c++ fmod is defined abnormally for negative numbers)
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE Number
pmod(const Number &value, const Number &modulus)
{
  return std::fmod(std::fmod(value, modulus) + modulus, modulus);
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
inline DEAL_II_ALWAYS_INLINE void
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
inline DEAL_II_ALWAYS_INLINE void
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
 * @brief Stabilization term as defined by
 * [Tezduyar](https://doi.org/10.1016/0045-7825(92)90141-6)
 *
 * For the element size, we use the diameter of a sphere with equivalent volume as the
 * cell. This is then divided by the degree.
 */
template <unsigned int dim, unsigned int degree, typename T>
inline DEAL_II_ALWAYS_INLINE T
stabilization_parameter(const T                         &timestep,
                        const T                         &element_volume,
                        const dealii::Tensor<1, dim, T> &velocity,
                        const T                         &kinematic_viscosity)
{
  using dealii::Utilities::fixed_power;
  using std::cbrt;
  using std::sqrt;

  constexpr auto degree_modifier  = T(1) / T(degree);
  constexpr auto inv_pi           = T(1) / T(M_PI);
  const auto     velocity_l2_norm = velocity.norm_square() + T(1e-12);

  auto element_size = T(0);
  if constexpr (dim == 1)
    {
      element_size = element_volume;
    }
  else if constexpr (dim == 2)
    {
      element_size = T(2) * sqrt(element_volume * inv_pi);
    }
  else if constexpr (dim == 3)
    {
      element_size = cbrt(T(6) * element_volume * inv_pi);
    }
  else
    {
      Assert(false, dealii::ExcMessage("Invalid dimension for the stabilization term."));
    }
  element_size *= degree_modifier;
  const auto element_size_2     = element_size * element_size;
  const auto element_size_2_inv = T(1) / element_size_2;
  const auto element_size_4_inv = element_size_2_inv * element_size_2_inv;

  const T term_time = fixed_power<2>(T(2) / timestep);
  const T term_adv  = T(4) * velocity_l2_norm * element_size_2_inv;
  const T term_diff = T(144) * fixed_power<2>(kinematic_viscosity) * element_size_4_inv;

  return T(1) / sqrt(term_time + term_adv + term_diff);
}

/**
 * @brief Stabilization term as defined by
 * [Tezduyar](https://doi.org/10.1016/0045-7825(92)90141-6)
 *
 * For the element size, we use the diameter of a sphere with equivalent volume as the
 * cell. This is then divided by the degree.
 *
 * @note This is for steady state problems
 */
template <unsigned int dim, unsigned int degree, typename T>
inline DEAL_II_ALWAYS_INLINE T
stabilization_parameter(const T                         &element_volume,
                        const dealii::Tensor<1, dim, T> &velocity,
                        const T                         &kinematic_viscosity)
{
  using dealii::Utilities::fixed_power;
  using std::cbrt;
  using std::sqrt;

  constexpr auto degree_modifier  = T(1) / T(degree);
  constexpr auto inv_pi           = T(1) / T(M_PI);
  const auto     velocity_l2_norm = velocity.norm_square() + T(1e-12);

  auto element_size = T(0);
  if constexpr (dim == 1)
    {
      element_size = element_volume;
    }
  else if constexpr (dim == 2)
    {
      element_size = T(2) * sqrt(element_volume * inv_pi);
    }
  else if constexpr (dim == 3)
    {
      element_size = cbrt(T(6) * element_volume * inv_pi);
    }
  else
    {
      Assert(false, dealii::ExcMessage("Invalid dimension for the stabilization term."));
    }
  element_size *= degree_modifier;
  const auto element_size_2     = element_size * element_size;
  const auto element_size_2_inv = T(1) / element_size_2;
  const auto element_size_4_inv = element_size_2_inv * element_size_2_inv;

  const T term_adv  = T(4) * velocity_l2_norm * element_size_2_inv;
  const T term_diff = T(144) * fixed_power<2>(kinematic_viscosity) * element_size_4_inv;

  return T(1) / sqrt(term_adv + term_diff);
}

/**
 * @brief Stabilization term as defined by
 * [Tezduyar](https://doi.org/10.1016/0045-7825(92)90141-6)
 *
 * For the element size, we use the diameter of a sphere with equivalent volume as the
 * cell. This is then divided by the degree.
 *
 * @note This is for pure advection transient problems
 */
template <unsigned int dim, unsigned int degree, typename T>
inline DEAL_II_ALWAYS_INLINE T
stabilization_parameter(const T                         &timestep,
                        const T                         &element_volume,
                        const dealii::Tensor<1, dim, T> &velocity)
{
  using dealii::Utilities::fixed_power;
  using std::cbrt;
  using std::sqrt;

  constexpr auto degree_modifier  = T(1) / T(degree);
  constexpr auto inv_pi           = T(1) / T(M_PI);
  const auto     velocity_l2_norm = velocity.norm_square() + T(1e-12);

  auto element_size = T(0);
  if constexpr (dim == 1)
    {
      element_size = element_volume;
    }
  else if constexpr (dim == 2)
    {
      element_size = T(2) * sqrt(element_volume * inv_pi);
    }
  else if constexpr (dim == 3)
    {
      element_size = cbrt(T(6) * element_volume * inv_pi);
    }
  else
    {
      Assert(false, dealii::ExcMessage("Invalid dimension for the stabilization term."));
    }
  element_size *= degree_modifier;
  const auto element_size_2     = element_size * element_size;
  const auto element_size_2_inv = T(1) / element_size_2;
  const auto element_size_4_inv = element_size_2_inv * element_size_2_inv;

  const T term_time = fixed_power<2>(T(2) / timestep);
  const T term_adv  = T(4) * velocity_l2_norm * element_size_2_inv;

  return T(1) / sqrt(term_time + term_adv);
}

/**
 * @brief Stabilization term as defined by
 * [Tezduyar](https://doi.org/10.1016/0045-7825(92)90141-6)
 *
 * For the element size, we use the diameter of a sphere with equivalent volume as the
 * cell. This is then divided by the degree.
 *
 * @note This is for pure advection steady-state problems
 */
template <unsigned int dim, unsigned int degree, typename T>
inline DEAL_II_ALWAYS_INLINE T
stabilization_parameter(const T                         &element_volume,
                        const dealii::Tensor<1, dim, T> &velocity)
{
  using dealii::Utilities::fixed_power;
  using std::cbrt;
  using std::sqrt;

  constexpr auto degree_modifier  = T(1) / T(degree);
  constexpr auto inv_pi           = T(1) / T(M_PI);
  const auto     velocity_l2_norm = velocity.norm_square() + T(1e-12);

  auto element_size = T(0);
  if constexpr (dim == 1)
    {
      element_size = element_volume;
    }
  else if constexpr (dim == 2)
    {
      element_size = T(2) * sqrt(element_volume * inv_pi);
    }
  else if constexpr (dim == 3)
    {
      element_size = cbrt(T(6) * element_volume * inv_pi);
    }
  else
    {
      Assert(false, dealii::ExcMessage("Invalid dimension for the stabilization term."));
    }
  element_size *= degree_modifier;
  const auto element_size_2     = element_size * element_size;
  const auto element_size_2_inv = T(1) / element_size_2;
  const auto element_size_4_inv = element_size_2_inv * element_size_2_inv;

  const T term_adv = T(4) * velocity_l2_norm * element_size_2_inv;
  // NOTE: I could simplify the math here, but I won't bother for now since I imagine the
  // use case is pretty limited.

  return T(1) / sqrt(term_adv);
}

/**
 * @brief Remove whitespace from strings
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
inline DEAL_II_ALWAYS_INLINE std::vector<number>
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
