// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <prismspf/core/type_enums.h>

#include <prismspf/utilities/logger.h>

#include <prismspf/config.h>

#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

namespace Mechanics
{
  /**
   * @brief Validate supported combinations (compile time).
   */
  template <unsigned int dim, StressState state>
  inline constexpr bool valid_stress_state =
    (dim == 1 && state == StressState::ThreeDimensional) ||
    (dim == 2 &&
     (state == StressState::PlaneStress || state == StressState::PlaneStrain)) ||
    (dim == 3 && state == StressState::ThreeDimensional);

  /**
   * @brief Voigt notation index range.
   * This is evaluated at compile time. The user must state the 2D assumption explicitly.
   *
   * get_voigt_size<1>();              // ThreeDimensional, returns 1
   * get_voigt_size<3>();              // ThreeDimensional, returns 6
   * get_voigt_size<2, StressState::PlaneStress>(); // returns 3
   * get_voigt_size<2, StressState::PlaneStrain>(); // returns 4
   * get_voigt_size<2>();              // this is invalid
   */
  template <unsigned int dim, StressState state = StressState::ThreeDimensional>
  constexpr unsigned int
  get_voigt_size()
  {
    static_assert(dim >= 1 && dim <= 3,
                  "Mechanics supports only dimensions 1, 2, and 3.");
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination");

    if constexpr (dim == 1)
      return 1;
    else if constexpr (dim == 2 && state == StressState::PlaneStress)
      return 3;
    else if constexpr (dim == 2 && state == StressState::PlaneStrain)
      return 4;
    else
      return 6;
  }

  /**
   * @brief Strain tensor to Voigt notation.
   * 1D, 2D Plane Stress, 3D.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T         = double>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE void
  strain_to_voigt(const dealii::Tensor<2, dim, T>                    &tensor,
                  dealii::Tensor<1, get_voigt_size<dim, state>(), T> &voigt)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    if constexpr (dim == 1)
      {
        voigt[0] = tensor[0][0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[0][1] + tensor[1][0];
      }
    else if constexpr (dim == 3)
      {
        // Voigt indexing (11, 22, 33, 23, 13, 12)
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[2][2];
        voigt[3] = tensor[1][2] + tensor[2][1];
        voigt[4] = tensor[0][2] + tensor[2][0];
        voigt[5] = tensor[0][1] + tensor[1][0];
      }
  }

  /**
   * @brief Strain tensor to Voigt notation.
   * 1D, 2D Plane Stress, 3D.
   * Overload: Return value
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T         = double>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<1, get_voigt_size<dim, state>(), T>
  strain_to_voigt(const dealii::Tensor<2, dim, T> &tensor)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    dealii::Tensor<1, get_voigt_size<dim, state>(), T> voigt;

    if constexpr (dim == 1)
      {
        voigt[0] = tensor[0][0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[0][1] + tensor[1][0];
      }
    else if constexpr (dim == 3)
      {
        // Voigt indexing (11, 22, 33, 23, 13, 12)
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[2][2];
        voigt[3] = tensor[1][2] + tensor[2][1];
        voigt[4] = tensor[0][2] + tensor[2][0];
        voigt[5] = tensor[0][1] + tensor[1][0];
      }

    return voigt;
  }

  /**
   * @brief Strain tensor to Voigt notation.
   * 1D, 2D Plane Stress, 3D.
   * Overload: Tensor input is a symmetric tensor.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T         = double>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE void
  strain_to_voigt(const dealii::SymmetricTensor<2, dim, T>           &tensor,
                  dealii::Tensor<1, get_voigt_size<dim, state>(), T> &voigt)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    if constexpr (dim == 1)
      {
        voigt[0] = tensor[0][0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = 2.0 * tensor[0][1];
      }
    else if constexpr (dim == 3)
      {
        // Voigt indexing (11, 22, 33, 23, 13, 12)
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[2][2];
        voigt[3] = 2.0 * tensor[1][2];
        voigt[4] = 2.0 * tensor[0][2];
        voigt[5] = 2.0 * tensor[0][1];
      }
  }

  /**
   * @brief Strain tensor to Voigt notation.
   * 1D, 2D Plane Stress, 3D.
   * Overload: Tensor input is a symmetric tensor.
   * Overload: Return value.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T         = double>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<1, get_voigt_size<dim, state>(), T>
  strain_to_voigt(const dealii::SymmetricTensor<2, dim, T> &tensor)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    dealii::Tensor<1, get_voigt_size<dim, state>(), T> voigt;

    if constexpr (dim == 1)
      {
        voigt[0] = tensor[0][0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = 2.0 * tensor[0][1];
      }
    else if constexpr (dim == 3)
      {
        // Voigt indexing (11, 22, 33, 23, 13, 12)
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[2][2];
        voigt[3] = 2.0 * tensor[1][2];
        voigt[4] = 2.0 * tensor[0][2];
        voigt[5] = 2.0 * tensor[0][1];
      }

    return voigt;
  }

  /**
   * @brief Strain tensor to Voigt notation.
   * Overload for 2D Plane Strain
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE void
  strain_to_voigt(const dealii::Tensor<2, dim, T> &tensor_inplane,
                  const T                         &component_zz,
                  dealii::Tensor<1, 4, T>         &voigt)
  {
    // Plane Strain
    voigt[0] = tensor_inplane[0][0];
    voigt[1] = tensor_inplane[1][1];
    voigt[2] = component_zz;
    voigt[3] = tensor_inplane[0][1] + tensor_inplane[1][0];
  }

  /**
   * @brief Strain tensor to Voigt notation.
   * Overload for 2D Plane Strain
   * Overload: Return value
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<1, 4, T>
  strain_to_voigt(const dealii::Tensor<2, dim, T> &tensor_inplane, const T &component_zz)
  {
    dealii::Tensor<1, 4, T> voigt;

    // Plane Strain
    voigt[0] = tensor_inplane[0][0];
    voigt[1] = tensor_inplane[1][1];
    voigt[2] = component_zz;
    voigt[3] = tensor_inplane[0][1] + tensor_inplane[1][0];

    return voigt;
  }

  /**
   * @brief Strain tensor to Voigt notation.
   * Overload for 2D Plane Strain
   * Overload: Tensor input is a symmetric tensor.
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE void
  strain_to_voigt(const dealii::SymmetricTensor<2, dim, T> &tensor_inplane,
                  const T                                  &component_zz,
                  dealii::Tensor<1, 4, T>                  &voigt)
  {
    // Plane Strain
    voigt[0] = tensor_inplane[0][0];
    voigt[1] = tensor_inplane[1][1];
    voigt[2] = component_zz;
    voigt[3] = 2.0 * tensor_inplane[0][1];
  }

  /**
   * @brief Strain tensor to Voigt notation.
   * Overload for 2D Plane Strain
   * Overload: Tensor input is a symmetric tensor.
   * Overload: Return value
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<1, 4, T>
  strain_to_voigt(const dealii::SymmetricTensor<2, dim, T> &tensor_inplane,
                  const T                                  &component_zz)
  {
    dealii::Tensor<1, 4, T> voigt;

    // Plane Strain
    voigt[0] = tensor_inplane[0][0];
    voigt[1] = tensor_inplane[1][1];
    voigt[2] = component_zz;
    voigt[3] = 2.0 * tensor_inplane[0][1];

    return voigt;
  }

  /**
   * @brief Voigt notation to Strain tensor.
   * 1D, 2D Plane Stress, 3D.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T         = double>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE void
  voigt_to_strain(const dealii::Tensor<1, get_voigt_size<dim, state>(), T> &voigt,
                  dealii::Tensor<2, dim, T>                                &tensor)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    if constexpr (dim == 1)
      {
        tensor[0][0] = voigt[0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[0][1] = tensor[1][0] = 0.5 * voigt[2];
      }
    else if constexpr (dim == 3)
      {
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[2][2] = voigt[2];
        tensor[1][2] = tensor[2][1] = 0.5 * voigt[3];
        tensor[0][2] = tensor[2][0] = 0.5 * voigt[4];
        tensor[0][1] = tensor[1][0] = 0.5 * voigt[5];
      }
  }

  /**
   * @brief Voigt notation to Strain tensor.
   * 1D, 2D Plane Stress, 3D.
   * Overload: Tensor output is a symmetric tensor.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T         = double>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE void
  voigt_to_strain(const dealii::Tensor<1, get_voigt_size<dim, state>(), T> &voigt,
                  dealii::SymmetricTensor<2, dim, T>                       &tensor)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    if constexpr (dim == 1)
      {
        tensor[0][0] = voigt[0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[0][1] = 0.5 * voigt[2];
      }
    else if constexpr (dim == 3)
      {
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[2][2] = voigt[2];
        tensor[1][2] = 0.5 * voigt[3];
        tensor[0][2] = 0.5 * voigt[4];
        tensor[0][1] = 0.5 * voigt[5];
      }
  }

  /**
   * @brief Voigt notation to Strain tensor.
   * 1D, 2D Plane Stress, 3D.
   * Overload: Return value, always return a symmetric tensor.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T         = double>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE dealii::SymmetricTensor<2, dim, T>
  voigt_to_strain(const dealii::Tensor<1, get_voigt_size<dim, state>(), T> &voigt)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    dealii::SymmetricTensor<2, dim, T> tensor;

    if constexpr (dim == 1)
      {
        tensor[0][0] = voigt[0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[0][1] = 0.5 * voigt[2];
      }
    else if constexpr (dim == 3)
      {
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[2][2] = voigt[2];
        tensor[1][2] = 0.5 * voigt[3];
        tensor[0][2] = 0.5 * voigt[4];
        tensor[0][1] = 0.5 * voigt[5];
      }

    return tensor;
  }

  /**
   * @brief Voigt notation to Strain tensor.
   * Overload for 2D Plane Strain
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE void
  voigt_to_strain(const dealii::Tensor<1, 4, T> &voigt,
                  dealii::Tensor<2, dim, T>     &tensor_inplane,
                  T                             &component_zz)
  {
    tensor_inplane[0][0] = voigt[0];
    tensor_inplane[1][1] = voigt[1];
    tensor_inplane[0][1] = tensor_inplane[1][0] = 0.5 * voigt[3];

    component_zz = voigt[2];
  }

  /**
   * @brief Voigt notation to Strain tensor.
   * Overload for 2D Plane Strain
   * Overload: Tensor output is a symmetric tensor.
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE void
  voigt_to_strain(const dealii::Tensor<1, 4, T>      &voigt,
                  dealii::SymmetricTensor<2, dim, T> &tensor_inplane,
                  T                                  &component_zz)
  {
    tensor_inplane[0][0] = voigt[0];
    tensor_inplane[1][1] = voigt[1];
    tensor_inplane[0][1] = 0.5 * voigt[3];

    component_zz = voigt[2];
  }

  /**
   * @brief Voigt notation to Strain tensor.
   * Overload for 2D Plane Strain
   * Overload: Return value, always return a symmetric tensor.
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE std::pair<dealii::SymmetricTensor<2, dim, T>, T>
                               voigt_to_strain(const dealii::Tensor<1, 4, T> &voigt)
  {
    dealii::SymmetricTensor<2, dim, T> tensor_inplane;

    tensor_inplane[0][0] = voigt[0];
    tensor_inplane[1][1] = voigt[1];
    tensor_inplane[0][1] = 0.5 * voigt[3];

    T component_zz = voigt[2];

    return {tensor_inplane, component_zz};
  }

  /**
   * @brief Stress tensor to Voigt notation.
   * 1D, 2D Plane Stress, 3D.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE void
  stress_to_voigt(const dealii::Tensor<2, dim, T>                    &tensor,
                  dealii::Tensor<1, get_voigt_size<dim, state>(), T> &voigt)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    if constexpr (dim == 1)
      {
        voigt[0] = tensor[0][0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = 0.5 * (tensor[0][1] + tensor[1][0]);
      }
    else if constexpr (dim == 3)
      {
        // Voigt indexing (11, 22, 33, 23, 13, 12)
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[2][2];
        voigt[3] = 0.5 * (tensor[1][2] + tensor[2][1]);
        voigt[4] = 0.5 * (tensor[0][2] + tensor[2][0]);
        voigt[5] = 0.5 * (tensor[0][1] + tensor[1][0]);
      }
  }

  /**
   * @brief Stress tensor to Voigt notation.
   * 1D, 2D Plane Stress, 3D.
   * Overload: Return value.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<1, get_voigt_size<dim, state>(), T>
  stress_to_voigt(const dealii::Tensor<2, dim, T> &tensor)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    dealii::Tensor<1, get_voigt_size<dim, state>(), T> voigt;

    if constexpr (dim == 1)
      {
        voigt[0] = tensor[0][0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = 0.5 * (tensor[0][1] + tensor[1][0]);
      }
    else if constexpr (dim == 3)
      {
        // Voigt indexing (11, 22, 33, 23, 13, 12)
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[2][2];
        voigt[3] = 0.5 * (tensor[1][2] + tensor[2][1]);
        voigt[4] = 0.5 * (tensor[0][2] + tensor[2][0]);
        voigt[5] = 0.5 * (tensor[0][1] + tensor[1][0]);
      }

    return voigt;
  }

  /**
   * @brief Stress tensor to Voigt notation.
   * 1D, 2D Plane Stress, 3D.
   * Overload: Tensor input is a symmetric tensor.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE void
  stress_to_voigt(const dealii::SymmetricTensor<2, dim, T>           &tensor,
                  dealii::Tensor<1, get_voigt_size<dim, state>(), T> &voigt)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    if constexpr (dim == 1)
      {
        voigt[0] = tensor[0][0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[0][1];
      }
    else if constexpr (dim == 3)
      {
        // Voigt indexing (11, 22, 33, 23, 13, 12)
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[2][2];
        voigt[3] = tensor[1][2];
        voigt[4] = tensor[0][2];
        voigt[5] = tensor[0][1];
      }
  }

  /**
   * @brief Stress tensor to Voigt notation.
   * 1D, 2D Plane Stress, 3D.
   * Overload: Tensor input is a symmetric tensor.
   * Overload: Return value.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<1, get_voigt_size<dim, state>(), T>
  stress_to_voigt(const dealii::SymmetricTensor<2, dim, T> &tensor)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    dealii::Tensor<1, get_voigt_size<dim, state>(), T> voigt;

    if constexpr (dim == 1)
      {
        voigt[0] = tensor[0][0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[0][1];
      }
    else if constexpr (dim == 3)
      {
        // Voigt indexing (11, 22, 33, 23, 13, 12)
        voigt[0] = tensor[0][0];
        voigt[1] = tensor[1][1];
        voigt[2] = tensor[2][2];
        voigt[3] = tensor[1][2];
        voigt[4] = tensor[0][2];
        voigt[5] = tensor[0][1];
      }

    return voigt;
  }

  /**
   * @brief Stress tensor to Voigt notation.
   * Overload for 2D Plane Strain
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE void
  stress_to_voigt(const dealii::Tensor<2, dim, T> &tensor_inplane,
                  const T                         &component_zz,
                  dealii::Tensor<1, 4, T>         &voigt)
  {
    // Plane Strain
    voigt[0] = tensor_inplane[0][0];
    voigt[1] = tensor_inplane[1][1];
    voigt[2] = component_zz;
    voigt[3] = tensor_inplane[0][1];
  }

  /**
   * @brief Stress tensor to Voigt notation.
   * Overload for 2D Plane Strain
   * Overload: Return value
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<1, 4, T>
  stress_to_voigt(const dealii::Tensor<2, dim, T> &tensor_inplane, const T &component_zz)
  {
    dealii::Tensor<1, 4, T> voigt;

    // Plane Strain
    voigt[0] = tensor_inplane[0][0];
    voigt[1] = tensor_inplane[1][1];
    voigt[2] = component_zz;
    voigt[3] = tensor_inplane[0][1];

    return voigt;
  }

  /**
   * @brief Stress tensor to Voigt notation.
   * Overload for 2D Plane Strain
   * Overload: Tensor input is a symmetric tensor.
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE void
  stress_to_voigt(const dealii::SymmetricTensor<2, dim, T> &tensor_inplane,
                  const T                                  &component_zz,
                  dealii::Tensor<1, 4, T>                  &voigt)
  {
    // Plane Strain
    voigt[0] = tensor_inplane[0][0];
    voigt[1] = tensor_inplane[1][1];
    voigt[2] = component_zz;
    voigt[3] = tensor_inplane[0][1];
  }

  /**
   * @brief Stress tensor to Voigt notation.
   * Overload for 2D Plane Strain
   * Overload: Tensor input is a symmetric tensor.
   * Overload: Return value
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<1, 4, T>
  stress_to_voigt(const dealii::SymmetricTensor<2, dim, T> &tensor_inplane,
                  const T                                  &component_zz)
  {
    dealii::Tensor<1, 4, T> voigt;

    // Plane Strain
    voigt[0] = tensor_inplane[0][0];
    voigt[1] = tensor_inplane[1][1];
    voigt[2] = component_zz;
    voigt[3] = tensor_inplane[0][1];

    return voigt;
  }

  /**
   * @brief Voigt notation to Stress tensor.
   * 1D, 2D Plane Stress, 3D.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T         = double>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE void
  voigt_to_stress(const dealii::Tensor<1, get_voigt_size<dim, state>(), T> &voigt,
                  dealii::Tensor<2, dim, T>                                &tensor)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    if constexpr (dim == 1)
      {
        tensor[0][0] = voigt[0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[0][1] = tensor[1][0] = voigt[2];
      }
    else if constexpr (dim == 3)
      {
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[2][2] = voigt[2];
        tensor[1][2] = tensor[2][1] = voigt[3];
        tensor[0][2] = tensor[2][0] = voigt[4];
        tensor[0][1] = tensor[1][0] = voigt[5];
      }
  }

  /**
   * @brief Voigt notation to Stress tensor.
   * 1D, 2D Plane Stress, 3D.
   * Overload: Tensor output is a symmetric tensor.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T         = double>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE void
  voigt_to_stress(const dealii::Tensor<1, get_voigt_size<dim, state>(), T> &voigt,
                  dealii::SymmetricTensor<2, dim, T>                       &tensor)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    if constexpr (dim == 1)
      {
        tensor[0][0] = voigt[0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[0][1] = voigt[2];
      }
    else if constexpr (dim == 3)
      {
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[2][2] = voigt[2];
        tensor[1][2] = voigt[3];
        tensor[0][2] = voigt[4];
        tensor[0][1] = voigt[5];
      }
  }

  /**
   * @brief Voigt notation to Stress tensor.
   * 1D, 2D Plane Stress, 3D.
   * Overload: Return value, always return a symmetric tensor.
   */
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T         = double>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE dealii::SymmetricTensor<2, dim, T>
  voigt_to_stress(const dealii::Tensor<1, get_voigt_size<dim, state>(), T> &voigt)
  {
    static_assert(valid_stress_state<dim, state>,
                  "Invalid dimension/StressState combination.");

    dealii::SymmetricTensor<2, dim, T> tensor;

    if constexpr (dim == 1)
      {
        tensor[0][0] = voigt[0];
      }
    else if constexpr (dim == 2)
      {
        // Plane Stress
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[0][1] = voigt[2];
      }
    else if constexpr (dim == 3)
      {
        tensor[0][0] = voigt[0];
        tensor[1][1] = voigt[1];
        tensor[2][2] = voigt[2];
        tensor[1][2] = voigt[3];
        tensor[0][2] = voigt[4];
        tensor[0][1] = voigt[5];
      }

    return tensor;
  }

  /**
   * @brief Voigt to Stress tensor.
   * Overload for 2D Plane Strain
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE void
  voigt_to_stress(const dealii::Tensor<1, 4, T> &voigt,
                  dealii::Tensor<2, dim, T>     &tensor_inplane,
                  T                             &component_zz)
  {
    tensor_inplane[0][0] = voigt[0];
    tensor_inplane[1][1] = voigt[1];
    tensor_inplane[0][1] = tensor_inplane[1][0] = voigt[3];

    component_zz = voigt[2];
  }

  /**
   * @brief Voigt to Stress tensor.
   * Overload for 2D Plane Strain
   * Overload: Tensor output is a symmetric tensor.
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE void
  voigt_to_stress(const dealii::Tensor<1, 4, T>      &voigt,
                  dealii::SymmetricTensor<2, dim, T> &tensor_inplane,
                  T                                  &component_zz)
  {
    tensor_inplane[0][0] = voigt[0];
    tensor_inplane[1][1] = voigt[1];
    tensor_inplane[0][1] = voigt[3];

    component_zz = voigt[2];
  }

  /**
   * @brief Voigt to Stress tensor.
   * Overload for 2D Plane Strain
   * Overload: Return value, always return a symmetric tensor.
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE std::pair<dealii::SymmetricTensor<2, dim, T>, T>
                               voigt_to_stress(const dealii::Tensor<1, 4, T> &voigt)
  {
    dealii::SymmetricTensor<2, dim, T> tensor_inplane;

    tensor_inplane[0][0] = voigt[0];
    tensor_inplane[1][1] = voigt[1];
    tensor_inplane[0][1] = voigt[3];

    T component_zz = voigt[2];

    return {tensor_inplane, component_zz};
  }

  /**
   * @brief Isotropic stiffness matrix.
   */
  // TODO: should we use DEAL_II_ALWAYS_INLINE
  template <unsigned int dim,
            StressState  state = StressState::ThreeDimensional,
            typename T         = double>
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<2, get_voigt_size<dim, state>(), T>
                               stiffness_isotropic(const T E, const T nu)
  {
    AssertThrow(E > 0,
                dealii::ExcMessage("Invalid isotropic elastic constants: "
                                   "Young's modulus E must be positive."));
    AssertThrow(nu > -1.0 && nu < 0.5,
                dealii::ExcMessage("Invalid isotropic elastic constants: "
                                   "Poisson's ratio must be in range -1 < nu < 0.5"));

    constexpr unsigned int           voigt_size = get_voigt_size<dim, state>();
    dealii::Tensor<2, voigt_size, T> stiffness;

    if constexpr (dim == 1)
      {
        stiffness[0][0] = E;
      }
    else if constexpr (dim == 2)
      {
        const T G = E / (2.0 * (1.0 + nu));

        if constexpr (state == StressState::PlaneStress)
          {
            // 11, 22, 12
            const T lambda  = (nu * E) / (1.0 - nu * nu);
            stiffness[0][0] = stiffness[1][1] = lambda + 2.0 * G;
            stiffness[0][1] = stiffness[1][0] = lambda;
            stiffness[2][2]                   = G;
          }
        else if constexpr (state == StressState::PlaneStrain)
          {
            // TODO: Warning for parameer close to incompressible

            // 11, 22, 33, 12
            const T lambda  = (nu * E) / ((1.0 + nu) * (1.0 - 2.0 * nu));
            stiffness[0][0] = stiffness[1][1] = stiffness[2][2] = lambda + 2.0 * G;
            stiffness[0][1] = stiffness[1][0] = lambda;
            stiffness[0][2] = stiffness[2][0] = lambda;
            stiffness[1][2] = stiffness[2][1] = lambda;
            stiffness[3][3]                   = G;
          }
        else
          {
            AssertThrow(false, dealii::ExcMessage("Invalid stress state type for 2D"));
          }
      }
    else if constexpr (dim == 3)
      {
        // TODO: Warning for parameer close to incompressible

        // 11, 22, 33, 23, 13, 12
        const T G      = E / (2.0 * (1.0 + nu));
        const T lambda = (nu * E) / ((1.0 + nu) * (1.0 - 2.0 * nu));

        stiffness[0][0] = stiffness[1][1] = stiffness[2][2] = lambda + 2.0 * G;
        stiffness[0][1] = stiffness[1][0] = lambda;
        stiffness[0][2] = stiffness[2][0] = lambda;
        stiffness[1][2] = stiffness[2][1] = lambda;

        stiffness[3][3] = G;
        stiffness[4][4] = G;
        stiffness[5][5] = G;
      }
    else
      {
        AssertThrow(false, dealii::ExcMessage("Unsupported dimension"));
      }

    return stiffness;
  }

  /**
   * @brief Orthotropic stiffness matrix: Overload for Plane Stress.
   */
  // TODO: should we use DEAL_II_ALWAYS_INLINE
  template <typename T>
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<2, 3, T>
  stiffness_orthotropic(const T E1, const T E2, const T nu12, const T G12)
  {
    AssertThrow(E1 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: E1 must be positive."));

    AssertThrow(E2 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: E2 must be positive."));

    AssertThrow(G12 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: G12 must be positive."));

    dealii::Tensor<2, 3, T> stiffness;

    const T nu21 = nu12 * (E2 / E1);

    AssertThrow(
      1.0 > nu12 * nu21,
      dealii::ExcMessage(
        "Invalid orthotropic elastic constants: 1 - nu12*nu21 must be positive."));

    const T delta = 1.0 - nu12 * nu21;

    AssertThrow(
      delta > 0.0,
      dealii::ExcMessage(
        "Invalid orthotropic elastic constants: the determinant must be positive."));

    constexpr T tolerance = std::numeric_limits<T>::epsilon() * 1e4;
    if (delta <= tolerance)
      {
        Logger::instance()
          << LogFormatter::warning(
               "WARNING: Orthotropic elastic constants are nearly singular.")
          << std::endl;
      }

    const T inv_delta = 1.0 / delta;

    stiffness[0][0] = E1 * inv_delta;
    stiffness[1][1] = E2 * inv_delta;
    stiffness[0][1] = stiffness[1][0] = (nu12 * E2) * inv_delta;
    stiffness[2][2]                   = G12;

    return stiffness;
  }

  /**
   * @brief Orthotropic stiffness matrix: Overload for Plane Strain.
   */
  // TODO: should we use DEAL_II_ALWAYS_INLINE
  template <typename T>
  inline dealii::Tensor<2, 4, T>
  stiffness_orthotropic(const T E1,
                        const T E2,
                        const T E3,
                        const T nu12,
                        const T nu13,
                        const T nu23,
                        const T G12)
  {
    AssertThrow(E1 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: E1 must be positive."));

    AssertThrow(E2 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: E2 must be positive."));

    AssertThrow(E3 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: E3 must be positive."));

    AssertThrow(G12 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: G12 must be positive."));

    dealii::Tensor<2, 4, T> stiffness;

    const T nu21 = nu12 * (E2 / E1);
    const T nu31 = nu13 * (E3 / E1);
    const T nu32 = nu23 * (E3 / E2);

    AssertThrow(
      1.0 > nu12 * nu21,
      dealii::ExcMessage(
        "Invalid orthotropic elastic constants: 1 - nu12*nu21 must be positive."));

    AssertThrow(
      1.0 > nu13 * nu31,
      dealii::ExcMessage(
        "Invalid orthotropic elastic constants: 1 - nu13*nu31 must be positive."));

    AssertThrow(
      1.0 > nu23 * nu32,
      dealii::ExcMessage(
        "Invalid orthotropic elastic constants: 1 - nu23*nu32 must be positive."));

    const T delta =
      1.0 - (nu12 * nu21) - (nu23 * nu32) - (nu13 * nu31) - (2.0 * nu12 * nu23 * nu31);

    AssertThrow(
      delta > 0.0,
      dealii::ExcMessage(
        "Invalid orthotropic elastic constants: the determinant must be positive."));

    // TODO: warning for nearly singular

    const T inv_delta = 1.0 / delta;

    stiffness[0][0] = E1 * (1.0 - nu23 * nu32) * inv_delta;
    stiffness[1][1] = E2 * (1.0 - nu13 * nu31) * inv_delta;
    stiffness[2][2] = E3 * (T(1) - nu12 * nu21) * inv_delta;
    stiffness[0][1] = stiffness[1][0] = E1 * (nu21 + nu31 * nu23) * inv_delta;
    stiffness[0][2] = stiffness[2][0] = E1 * (nu31 + nu21 * nu32) * inv_delta;
    stiffness[1][2] = stiffness[2][1] = E2 * (nu32 + nu12 * nu31) * inv_delta;
    stiffness[3][3]                   = G12;

    return stiffness;
  }

  /**
   * @brief Orthotropic stiffness matrix: Overload for 3D.
   */
  // TODO: should we use DEAL_II_ALWAYS_INLINE
  template <typename T>
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<2, 6, T>
                               stiffness_orthotropic(const T E1,
                                                     const T E2,
                                                     const T E3,
                                                     const T nu12,
                                                     const T nu13,
                                                     const T nu23,
                                                     const T G12,
                                                     const T G13,
                                                     const T G23)
  {
    AssertThrow(E1 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: E1 must be positive."));

    AssertThrow(E2 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: E2 must be positive."));

    AssertThrow(E3 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: E3 must be positive."));

    AssertThrow(G12 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: G12 must be positive."));

    AssertThrow(G13 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: G13 must be positive."));

    AssertThrow(G23 > 0.0,
                dealii::ExcMessage(
                  "Invalid orthotropic elastic constants: G23 must be positive."));

    dealii::Tensor<2, 6, T> stiffness;

    const T nu21 = nu12 * (E2 / E1);
    const T nu31 = nu13 * (E3 / E1);
    const T nu32 = nu23 * (E3 / E2);

    AssertThrow(
      1.0 > nu12 * nu21,
      dealii::ExcMessage(
        "Invalid orthotropic elastic constants: 1 - nu12*nu21 must be positive."));

    AssertThrow(
      1.0 > nu13 * nu31,
      dealii::ExcMessage(
        "Invalid orthotropic elastic constants: 1 - nu13*nu31 must be positive."));

    AssertThrow(
      1.0 > nu23 * nu32,
      dealii::ExcMessage(
        "Invalid orthotropic elastic constants: 1 - nu23*nu32 must be positive."));

    const T delta =
      1.0 - (nu12 * nu21) - (nu23 * nu32) - (nu13 * nu31) - (2.0 * nu12 * nu23 * nu31);

    AssertThrow(
      delta > 0.0,
      dealii::ExcMessage(
        "Invalid orthotropic elastic constants: the determinant must be positive."));

    // TODO: warning for nearly singular

    const T inv_delta = 1.0 / delta;

    stiffness[0][0] = E1 * (1.0 - nu23 * nu32) * inv_delta;
    stiffness[1][1] = E2 * (1.0 - nu13 * nu31) * inv_delta;
    stiffness[2][2] = E3 * (1.0 - nu12 * nu21) * inv_delta;

    stiffness[0][1] = stiffness[1][0] = E1 * (nu21 + nu31 * nu23) * inv_delta;
    stiffness[0][2] = stiffness[2][0] = E1 * (nu31 + nu21 * nu32) * inv_delta;
    stiffness[1][2] = stiffness[2][1] = E2 * (nu32 + nu12 * nu31) * inv_delta;

    stiffness[3][3] = G23;
    stiffness[4][4] = G13;
    stiffness[5][5] = G12;

    return stiffness;
  }

  /**
   * @brief Extract 4x4 plane strain stiffness from 6x6 3D stiffness.
   * @note This restricted form assumes the xz and yz elastic shear components are zero,
   * which may not work for general anisotropy and eigenstrain.
   */
  template <typename T>
  inline DEAL_II_ALWAYS_INLINE dealii::Tensor<2, 4, T>
  extract_plane_strain_stiffness(const dealii::Tensor<2, 6, T> &stiffness_3d)
  {
    dealii::Tensor<2, 4, T> stiffness;

    constexpr std::array<unsigned int, 4> index = {0, 1, 2, 5};

    for (unsigned int i = 0; i < 4; ++i)
      for (unsigned int j = 0; j < 4; ++j)
        stiffness[i][j] = stiffness_3d[index[i]][index[j]];

    return stiffness;
  }

  /**
   * @brief Compute the stress with a given displacement and elasticity tensor. This
   * assumes that the provided parameters are in Voigt notation.
   * @note: Strain input is the elastic strain.
   */
  template <unsigned int dim, StressState state, typename T>
  inline DEAL_II_ALWAYS_INLINE void
  compute_stress(
    const dealii::Tensor<2, get_voigt_size<dim, state>(), T> &elasticity_tensor,
    const dealii::Tensor<1, get_voigt_size<dim, state>(), T> &strain,
    dealii::Tensor<1, get_voigt_size<dim, state>(), T>       &stress)
  {
    stress = elasticity_tensor * strain;
  }

  /**
   * @brief Compute the stress with a elasticity tensor (Voigt notation) and stress &
   * strain tensors.
   * 1D, 2D Plane Stress, 3D.
   *
   * @note This function internally converts to Voigt notation.
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state != StressState::PlaneStrain || dim != 2)
  inline DEAL_II_ALWAYS_INLINE void
  compute_stress(
    const dealii::Tensor<2, get_voigt_size<dim, state>(), T> &elasticity_tensor,
    const dealii::Tensor<2, dim, T>                          &strain,
    dealii::Tensor<2, dim, T>                                &stress)
  {
    dealii::Tensor<1, get_voigt_size<dim, state>(), T> sigma;
    dealii::Tensor<1, get_voigt_size<dim, state>(), T> epsilon;

    strain_to_voigt<dim, state, T>(strain, epsilon);
    compute_stress<dim, state, T>(elasticity_tensor, epsilon, sigma);
    voigt_to_stress<dim, state, T>(sigma, stress);
  }

  /**
   * @brief Compute the stress with a elasticity tensor (Voigt notation) and stress &
   * strain tensors.
   * Overload for 2D Plane Strain.
   * Input the in-plane strain and the out-of-plane component.
   * Return in-plane stress and out-of-plane stress component.
   *
   * @note This function internally converts to Voigt notation.
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE void
  compute_stress(const dealii::Tensor<2, 4, T>   &elasticity_tensor,
                 const dealii::Tensor<2, dim, T> &strain,
                 const T                         &strain_zz,
                 dealii::Tensor<2, dim, T>       &stress,
                 T                               &stress_zz)
  {
    dealii::Tensor<1, 4, T> sigma;
    dealii::Tensor<1, 4, T> epsilon;

    strain_to_voigt<dim, state, T>(strain, strain_zz, epsilon);
    compute_stress<dim, state, T>(elasticity_tensor, epsilon, sigma);
    voigt_to_stress<dim, state, T>(sigma, stress, stress_zz);
  }

  /**
   * @brief Compute the stress with a elasticity tensor and stress & strain tensors.
   * Overload for 2D Plane Strain.
   * Input the in-plane strain and the out-of-plane component.
   * Return stress in Voigt notation.
   *
   * @note This function internally converts to Voigt notation.
   */
  template <unsigned int dim, StressState state, typename T>
  requires(state == StressState::PlaneStrain && dim == 2)
  inline DEAL_II_ALWAYS_INLINE void
  compute_stress(const dealii::Tensor<2, 4, T>   &elasticity_tensor,
                 const dealii::Tensor<2, dim, T> &strain,
                 const T                         &strain_zz,
                 dealii::Tensor<1, 4, T>         &stress)
  {
    dealii::Tensor<1, 4, T> epsilon;

    strain_to_voigt<dim, state, T>(strain, strain_zz, epsilon);
    compute_stress<dim, state, T>(elasticity_tensor, epsilon, stress);
  }

  /* TODO: out-of-plane strain epsilon_zz under plane stress may be needed.*/

  /**
   * -------------------------------------------------------------
   * The following are provided just for backward compatibility.
   * It is preferred to use the functions above and provide StressState when calling
   * compute_stress.
   * -------------------------------------------------------------
   */

  /**
   * @brief Voigt notation index range.
   * This is used for backward compatibility.
   */
  template <unsigned int dim>
  [[deprecated("Use get_voigt_size<dim, StressState>() instead.")]] constexpr unsigned int
    voigt_tensor_size = (2 * dim) - 1 + (dim / 3);

  /**
   * @brief Compute the stress with a given displacement and elasticity tensor. This
   * assumes that the provided parameters are in Voigt notation.
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
   * @brief Compute the stress with a given displacement and elasticity tensor.
   *
   * @note This function internally converts to Voigt notation.
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

        // In Voigt notation: epsilon are engineering shear strains
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

        // In Voigt notation: epsilon are engineering shear strains
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

        stress[xx_dir][xx_dir] =
          elasticity_tensor[xx_dir][xx_dir] * strain[xx_dir][xx_dir];
      }
  }

} // namespace Mechanics

PRISMS_PF_END_NAMESPACE
