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

// NOLINTBEGIN(readability-identifier-naming,readability-magic-numbers,cppcoreguidelines-avoid-magic-numbers,readability-identifier-length)
PRISMS_PF_BEGIN_NAMESPACE
// Tolerance for checking singularity in elastic stiffness
template <typename T>
constexpr T tolerance = std::numeric_limits<T>::epsilon() * 1e4;

/**
 * @brief Class containing solid mechanics utility functions for a templated dimension.
 * 2D assumes plane stress conditions. For plane strain conditions, use the PlaneStrain
 * class instead.
 */
template <unsigned int dim>
struct Mechanics
{
  static_assert(dim >= 1 && dim <= 3, "Mechanics supports only dimensions 1, 2, and 3.");

  /**
   * @brief Voigt notation index range.
   * This is evaluated at compile time.
   *
   * dim = 1; // returns 1
   * dim = 2; // returns 3
   * dim = 3; // returns 6
   */
  static constexpr unsigned int voigt_size = (dim * dim - dim) / 2 + dim;

  template <typename T = double>
  using VoigtVector = dealii::Tensor<1, voigt_size, T>;
  template <typename T = double>
  using VoigtMatrix = dealii::Tensor<2, voigt_size, T>;
  template <typename T = double>
  using MechTensor = dealii::Tensor<2, dim, T>;
  template <typename T = double>
  using SymMechTensor = dealii::SymmetricTensor<2, dim, T>;

  /**
   * @brief Strain tensor to Voigt notation.
   * 1D, 2D, 3D.
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE void
  strain_to_voigt(const MechTensor<T> &tensor, VoigtVector<T> &voigt)
  {
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
   * 1D, 2D, 3D.
   * Overload: Return value
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE VoigtVector<T>
                                      strain_to_voigt(const MechTensor<T> &tensor)
  {
    VoigtVector<T> voigt;

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
   * 1D, 2D, 3D.
   * Overload: Tensor input is a symmetric tensor.
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE void
  strain_to_voigt(const SymMechTensor<T> &tensor, VoigtVector<T> &voigt)
  {
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
   * 1D, 2D, 3D.
   * Overload: Tensor input is a symmetric tensor.
   * Overload: Return value.
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE VoigtVector<T>
                                      strain_to_voigt(const SymMechTensor<T> &tensor)
  {
    VoigtVector<T> voigt;

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
   * @brief Voigt notation to Strain tensor.
   * 1D, 2D, 3D.
   */
  template <typename T = double>

  static inline DEAL_II_ALWAYS_INLINE void
  voigt_to_strain(const VoigtVector<T> &voigt, MechTensor<T> &tensor)
  {
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
   * 1D, 2D, 3D.
   * Overload: Tensor output is a symmetric tensor.
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE void
  voigt_to_strain(const VoigtVector<T> &voigt, SymMechTensor<T> &tensor)
  {
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
   * 1D, 2D, 3D.
   * Overload: Return value, always return a symmetric tensor.
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE SymMechTensor<T>
                                      voigt_to_strain(const VoigtVector<T> &voigt)
  {
    SymMechTensor<T> tensor;

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
   * @brief Stress tensor to Voigt notation.
   * 1D, 2D, 3D.
   */
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE void
  stress_to_voigt(const MechTensor<T> &tensor, VoigtVector<T> &voigt)
  {
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
   * 1D, 2D, 3D.
   * Overload: Return value.
   */
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE VoigtVector<T>
                                      stress_to_voigt(const MechTensor<T> &tensor)
  {
    VoigtVector<T> voigt;

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
   * 1D, 2D, 3D.
   * Overload: Tensor input is a symmetric tensor.
   */
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE void
  stress_to_voigt(const SymMechTensor<T> &tensor, VoigtVector<T> &voigt)
  {
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
   * 1D, 2D, 3D.
   * Overload: Tensor input is a symmetric tensor.
   * Overload: Return value.
   */
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE VoigtVector<T>
                                      stress_to_voigt(const SymMechTensor<T> &tensor)
  {
    VoigtVector<T> voigt;

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
   * @brief Voigt notation to Stress tensor.
   * 1D, 2D, 3D.
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE void
  voigt_to_stress(const VoigtVector<T> &voigt, MechTensor<T> &tensor)
  {
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
   * 1D, 2D, 3D.
   * Overload: Tensor output is a symmetric tensor.
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE void
  voigt_to_stress(const VoigtVector<T> &voigt, SymMechTensor<T> &tensor)
  {
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
   * 1D, 2D, 3D.
   * Overload: Return value, always return a symmetric tensor.
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE SymMechTensor<T>
                                      voigt_to_stress(const VoigtVector<T> &voigt)
  {
    SymMechTensor<T> tensor;

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
   * @brief Isotropic stiffness matrix.
   * 1D, 2D, 3D.
   */
  // TODO: should we use DEAL_II_ALWAYS_INLINE
  template <typename T = double>
  static inline VoigtMatrix<T>
  stiffness_isotropic(const T E, const T nu)
  {
    Assert(E > T(0.0),
           dealii::ExcMessage("Invalid isotropic elastic constants: "
                              "Young's modulus E must be positive."));
    Assert(nu > T(-1.0) && nu < T(0.5),
           dealii::ExcMessage("Invalid isotropic elastic constants: "
                              "Poisson's ratio must be in range -1 < nu < 0.5"));

    VoigtMatrix<T> stiffness;

    if constexpr (dim == 1)
      {
        stiffness[0][0] = E;
      }
    else if constexpr (dim == 2)
      {
        const T G = E / (T(2.0) * (T(1.0) + nu));

        // 11, 22, 12
        const T lambda  = (nu * E) / (T(1.0) - nu * nu);
        stiffness[0][0] = stiffness[1][1] = lambda + T(2.0) * G;
        stiffness[0][1] = stiffness[1][0] = lambda;
        stiffness[2][2]                   = G;
      }
    else if constexpr (dim == 3)
      {
#ifdef DEBUG
        // Warning for parameer close to incompressible
        if (std::fabs(T(1.0) - T(2.0) * nu) <= tolerance<T>)
          {
            Logger::instance()
              << LogFormatter::warning(
                   "WARNING: Isotropic elastic constants are nearly singular.")
              << std::endl;
          }
#endif

        // 11, 22, 33, 23, 13, 12
        const T G      = E / (T(2.0) * (T(1.0) + nu));
        const T lambda = (nu * E) / ((T(1.0) + nu) * (T(1.0) - T(2.0) * nu));

        stiffness[0][0] = stiffness[1][1] = stiffness[2][2] = lambda + T(2.0) * G;
        stiffness[0][1] = stiffness[1][0] = lambda;
        stiffness[0][2] = stiffness[2][0] = lambda;
        stiffness[1][2] = stiffness[2][1] = lambda;

        stiffness[3][3] = G;
        stiffness[4][4] = G;
        stiffness[5][5] = G;
      }

    return stiffness;
  }

  /**
   * @brief Orthotropic stiffness matrix.
   * 2D (Plane Stress)
   */
  // TODO: should we use DEAL_II_ALWAYS_INLINE
  template <typename T>
  requires(dim == 2)
  static inline VoigtMatrix<T>
  stiffness_orthotropic(const T E1, const T E2, const T nu12, const T G12)
  {
    Assert(E1 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: E1 must be positive."));

    Assert(E2 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: E2 must be positive."));

    Assert(G12 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: G12 must be positive."));

    VoigtMatrix<T> stiffness;

    const T nu21 = nu12 * (E2 / E1);

    Assert(T(1.0) > nu12 * nu21,
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: 1 - nu12*nu21 must be positive."));

    const T delta = T(1.0) - nu12 * nu21;

    Assert(delta > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: the determinant must be positive."));
#ifdef DEBUG
    if (delta <= tolerance<T>)
      {
        Logger::instance()
          << LogFormatter::warning(
               "WARNING: Orthotropic elastic constants are nearly singular.")
          << std::endl;
      }
#endif

    const T inv_delta = T(1.0) / delta;

    stiffness[0][0] = E1 * inv_delta;
    stiffness[1][1] = E2 * inv_delta;
    stiffness[0][1] = stiffness[1][0] = (nu12 * E2) * inv_delta;
    stiffness[2][2]                   = G12;

    return stiffness;
  }

  /**
   * @brief Orthotropic stiffness matrix.
   * 3D
   */
  // TODO: should we use DEAL_II_ALWAYS_INLINE
  template <typename T>
  requires(dim == 3)
  static inline VoigtMatrix<T>
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
    Assert(E1 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: E1 must be positive."));

    Assert(E2 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: E2 must be positive."));

    Assert(E3 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: E3 must be positive."));

    Assert(G12 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: G12 must be positive."));

    Assert(G13 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: G13 must be positive."));

    Assert(G23 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: G23 must be positive."));

    VoigtMatrix<T> stiffness;

    const T nu21 = nu12 * (E2 / E1);
    const T nu31 = nu13 * (E3 / E1);
    const T nu32 = nu23 * (E3 / E2);

    Assert(T(1.0) > nu12 * nu21,
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: 1 - nu12*nu21 must be positive."));

    Assert(T(1.0) > nu13 * nu31,
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: 1 - nu13*nu31 must be positive."));

    Assert(T(1.0) > nu23 * nu32,
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: 1 - nu23*nu32 must be positive."));

    const T delta = T(1.0) - (nu12 * nu21) - (nu23 * nu32) - (nu13 * nu31) -
                    (T(2.0) * nu12 * nu23 * nu31);

    Assert(delta > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: the determinant must be positive."));
#ifdef DEBUG
    // Warning for nearly singular
    if (delta <= tolerance<T>)
      {
        Logger::instance()
          << LogFormatter::warning(
               "WARNING: Orthotropic elastic constants are nearly singular.")
          << std::endl;
      }
#endif

    const T inv_delta = T(1.0) / delta;

    stiffness[0][0] = E1 * (T(1.0) - nu23 * nu32) * inv_delta;
    stiffness[1][1] = E2 * (T(1.0) - nu13 * nu31) * inv_delta;
    stiffness[2][2] = E3 * (T(1.0) - nu12 * nu21) * inv_delta;

    stiffness[0][1] = stiffness[1][0] = E1 * (nu21 + nu31 * nu23) * inv_delta;
    stiffness[0][2] = stiffness[2][0] = E1 * (nu31 + nu21 * nu32) * inv_delta;
    stiffness[1][2] = stiffness[2][1] = E2 * (nu32 + nu12 * nu31) * inv_delta;

    stiffness[3][3] = G23;
    stiffness[4][4] = G13;
    stiffness[5][5] = G12;

    return stiffness;
  }

  /**
   * @brief Compute the stress with a given displacement and elasticity tensor. This
   * assumes that the provided parameters are in Voigt notation.
   * @note: Strain input is the elastic strain.
   * 1D, 2D, 3D.
   */
  template <typename Tstiff, typename Tstrain>
  static inline DEAL_II_ALWAYS_INLINE auto // VoigtVector
  compute_stress(const VoigtMatrix<Tstiff>  &elasticity_tensor,
                 const VoigtVector<Tstrain> &strain)
  {
    return elasticity_tensor * strain;
  }

  /**
   * @brief Compute the stress with a given displacement and elasticity tensor. This
   * assumes that the provided parameters are in Voigt notation.
   * @note: Strain input is the elastic strain.
   * 1D, 2D, 3D.
   */
  template <typename Tstiff, typename Tstrain, typename Tstress>
  static inline DEAL_II_ALWAYS_INLINE void
  compute_stress(const VoigtMatrix<Tstiff>  &elasticity_tensor,
                 const VoigtVector<Tstrain> &strain,
                 VoigtVector<Tstress>       &stress)
  {
    stress = elasticity_tensor * strain;
  }

  /**
   * @brief Compute the stress with a elasticity tensor (Voigt notation) and stress &
   * strain tensors.
   * 2D, 3D.
   *
   * @note This function internally converts to Voigt notation.
   */
  template <typename Tstiff, typename Tstrain>
  requires(dim != 1)
  static inline DEAL_II_ALWAYS_INLINE auto // MechTensor
  compute_stress(const VoigtMatrix<Tstiff> &elasticity_tensor,
                 const MechTensor<Tstrain> &strain)
  {
    return voigt_to_stress(compute_stress(elasticity_tensor, strain_to_voigt(strain)));
  }

  /**
   * @brief Compute the stress with a elasticity tensor (Voigt notation) and stress &
   * strain tensors.
   * 2D, 3D.
   *
   * @note This function internally converts to Voigt notation.
   */
  template <typename Tstiff, typename Tstrain, typename Tstress>
  requires(dim != 1)
  static inline DEAL_II_ALWAYS_INLINE void
  compute_stress(const VoigtMatrix<Tstiff> &elasticity_tensor,
                 const MechTensor<Tstrain> &strain,
                 MechTensor<Tstress>       &stress)
  {
    stress = compute_stress(elasticity_tensor, strain);
  }

  /* TODO: out-of-plane strain epsilon_zz under plane stress may be needed. */

  /**
   * @brief Strain energy (Inputs are in Voigt notation).
   * 1D, 2D, 3D.
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE T
  strain_energy(const VoigtVector<T> &stress, const VoigtVector<T> &strain_e)
  {
    return T(0.5) * stress * strain_e;
  }

  /**
   * @brief von Mises stress (Input is in Voigt notation).
   * 1D, 2D, 3D.
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE T
  stress_mises(const VoigtVector<T> &stress)
  {
    T stress_m2 = T(0);

    if constexpr (dim == 1)
      {
        stress_m2 = stress[0] * stress[0];
      }
    else if constexpr (dim == 2)
      {
        const T &sigma_xx = stress[0];
        const T &sigma_yy = stress[1];
        const T &sigma_xy = stress[2];

        stress_m2 = sigma_xx * sigma_xx - sigma_xx * sigma_yy + sigma_yy * sigma_yy +
                    T(3) * sigma_xy * sigma_xy;
      }
    else if constexpr (dim == 3)
      {
        const T &sigma_xx = stress[0];
        const T &sigma_yy = stress[1];
        const T &sigma_zz = stress[2];
        const T &sigma_yz = stress[3];
        const T &sigma_xz = stress[4];
        const T &sigma_xy = stress[5];

        const T d_xy = sigma_xx - sigma_yy;
        const T d_yz = sigma_yy - sigma_zz;
        const T d_zx = sigma_zz - sigma_xx;

        stress_m2 =
          T(0.5) * (d_xy * d_xy + d_yz * d_yz + d_zx * d_zx) +
          T(3) * (sigma_xy * sigma_xy + sigma_yz * sigma_yz + sigma_xz * sigma_xz);
      }
    return std::sqrt(stress_m2);
  }

  /**
   * @brief Principal stress (Input is in Voigt notation).
   * @note For plane strain, returns the principal stresses of the in-plane 2x2 stress
   * tensor.
   * 1D, 2D.
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE dealii::Tensor<1, dim, T>
                                      stress_principal(const VoigtVector<T> &stress)
  {
    static_assert(
      dim == 1 || dim == 2,
      "stress_principal() supports only 1D, 2D plane stress, and 2D plane strain.");

    dealii::Tensor<1, dim, T> stress_p {};

    if constexpr (dim == 1)
      {
        stress_p[0] = stress[0];
      }
    else if constexpr (dim == 2)
      {
        constexpr unsigned int idx = 2;

        const T avg = (stress[0] + stress[1]) * T(0.5);
        const T dif = (stress[0] - stress[1]) * T(0.5);
        const T rad = std::sqrt(dif * dif + stress[idx] * stress[idx]);

        stress_p[0] = avg + rad;
        stress_p[1] = avg - rad;
      }
    else if constexpr (dim == 3)
      {
        // TODO: 3D
      }
    return stress_p;
  }

}; // class Mechanics

using PlaneStress = Mechanics<2>;

template <typename T = double>
using PlaneStressVector = PlaneStress::VoigtVector<T>;

template <typename T = double>
using PlaneStressStiffness = PlaneStress::VoigtMatrix<T>;

template <typename T = double>
using ThreeDimensionalVector = Mechanics<3>::VoigtVector<T>;

template <typename T = double>
using ThreeDimensionalStiffness = Mechanics<3>::VoigtMatrix<T>;

struct PlaneStrain
{
  /**
   * NOTE: The reduced 4x4 stiffness for plane strain requires stiffness like
   * [ .  .  .  0  0  . ]
   * [ .  .  .  0  0  . ]
   * [ .  .  .  0  0  . ]
   * [ 0  0  0  .  .  0 ]
   * [ 0  0  0  .  .  0 ]
   * [ .  .  .  0  0  . ]
   */
  static constexpr unsigned int dim        = 2;
  static constexpr unsigned int voigt_size = 4;
  template <typename T = double>
  using VoigtVector = dealii::Tensor<1, voigt_size, T>;
  template <typename T = double>
  using VoigtMatrix = dealii::Tensor<2, voigt_size, T>;
  template <typename T = double>
  using MechTensor = dealii::Tensor<2, dim, T>;
  template <typename T = double>
  using SymMechTensor = dealii::SymmetricTensor<2, dim, T>;

  /**
   * @brief Strain tensor to Voigt notation.
   * Overload for 2D Plane Strain
   */
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE void
  strain_to_voigt(const MechTensor<T> &tensor_inplane,
                  const T             &component_zz,
                  VoigtVector<T>      &voigt)
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
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE VoigtVector<T>
  strain_to_voigt(const MechTensor<T> &tensor_inplane, const T &component_zz)
  {
    VoigtVector<T> voigt;

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
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE void
  strain_to_voigt(const SymMechTensor<T> &tensor_inplane,
                  const T                &component_zz,
                  VoigtVector<T>         &voigt)
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
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE VoigtVector<T>
  strain_to_voigt(const SymMechTensor<T> &tensor_inplane, const T &component_zz)
  {
    VoigtVector<T> voigt;

    // Plane Strain
    voigt[0] = tensor_inplane[0][0];
    voigt[1] = tensor_inplane[1][1];
    voigt[2] = component_zz;
    voigt[3] = 2.0 * tensor_inplane[0][1];

    return voigt;
  }

  // -------------------------------------------------------------------------------

  /**
   * @brief Voigt notation to Strain tensor.
   * Overload for 2D Plane Strain
   */
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE void
  voigt_to_strain(const VoigtVector<T> &voigt,
                  MechTensor<T>        &tensor_inplane,
                  T                    &component_zz)
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
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE void
  voigt_to_strain(const VoigtVector<T> &voigt,
                  SymMechTensor<T>     &tensor_inplane,
                  T                    &component_zz)
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
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE std::pair<SymMechTensor<T>, T>
                                      voigt_to_strain(const VoigtVector<T> &voigt)
  {
    SymMechTensor<T> tensor_inplane;

    tensor_inplane[0][0] = voigt[0];
    tensor_inplane[1][1] = voigt[1];
    tensor_inplane[0][1] = 0.5 * voigt[3];

    T component_zz = voigt[2];

    return {tensor_inplane, component_zz};
  }

  // ------------------------------------------------------------------------

  /**
   * @brief Stress tensor to Voigt notation.
   * Overload for 2D Plane Strain
   */
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE void
  stress_to_voigt(const MechTensor<T> &tensor_inplane,
                  const T             &component_zz,
                  VoigtVector<T>      &voigt)
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
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE VoigtVector<T>
  stress_to_voigt(const MechTensor<T> &tensor_inplane, const T &component_zz)
  {
    VoigtVector<T> voigt;

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
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE void
  stress_to_voigt(const SymMechTensor<T> &tensor_inplane,
                  const T                &component_zz,
                  VoigtVector<T>         &voigt)
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
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE VoigtVector<T>
  stress_to_voigt(const SymMechTensor<T> &tensor_inplane, const T &component_zz)
  {
    VoigtVector<T> voigt;

    // Plane Strain
    voigt[0] = tensor_inplane[0][0];
    voigt[1] = tensor_inplane[1][1];
    voigt[2] = component_zz;
    voigt[3] = tensor_inplane[0][1];

    return voigt;
  }

  // ----------------------------------------------------------------------------

  /**
   * @brief Voigt to Stress tensor.
   * Overload for 2D Plane Strain
   */
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE void
  voigt_to_stress(const VoigtVector<T> &voigt,
                  MechTensor<T>        &tensor_inplane,
                  T                    &component_zz)
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
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE void
  voigt_to_stress(const VoigtVector<T> &voigt,
                  SymMechTensor<T>     &tensor_inplane,
                  T                    &component_zz)
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
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE std::pair<SymMechTensor<T>, T>
                                      voigt_to_stress(const VoigtVector<T> &voigt)
  {
    SymMechTensor<T> tensor_inplane;

    tensor_inplane[0][0] = voigt[0];
    tensor_inplane[1][1] = voigt[1];
    tensor_inplane[0][1] = voigt[3];

    T component_zz = voigt[2];

    return {tensor_inplane, component_zz};
  }

  //---------------------------------------------------------------

  /**
   * @brief Compute the stress with a given displacement and elasticity tensor. This
   * assumes that the provided parameters are in Voigt notation.
   * @note: Strain input is the elastic strain.
   */
  template <typename Tstiff, typename Tstrain>
  static inline DEAL_II_ALWAYS_INLINE auto // VoigtVector
  compute_stress(const VoigtMatrix<Tstiff>  &elasticity_tensor,
                 const VoigtVector<Tstrain> &strain)
  {
    return elasticity_tensor * strain;
  }

  /**
   * @brief Compute the stress with a given displacement and elasticity tensor. This
   * assumes that the provided parameters are in Voigt notation.
   * @note: Strain input is the elastic strain.
   */
  template <typename Tstiff, typename Tstrain, typename Tstress>
  static inline DEAL_II_ALWAYS_INLINE void
  compute_stress(const VoigtMatrix<Tstiff>  &elasticity_tensor,
                 const VoigtVector<Tstrain> &strain,
                 VoigtVector<Tstress>       &stress)
  {
    stress = elasticity_tensor * strain;
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
  template <typename Tstiff, typename Tstrain, typename Tstress>
  static inline DEAL_II_ALWAYS_INLINE void
  compute_stress(const VoigtMatrix<Tstiff> &elasticity_tensor,
                 const MechTensor<Tstrain> &strain,
                 const Tstrain             &strain_zz,
                 MechTensor<Tstress>       &stress,
                 Tstress                   &stress_zz)
  {
    VoigtVector<Tstrain> epsilon = strain_to_voigt(strain, strain_zz);
    VoigtVector<Tstress> sigma   = compute_stress(elasticity_tensor, epsilon, sigma);
    voigt_to_stress(sigma, stress, stress_zz);
  }

  /**
   * @brief Compute the stress with a elasticity tensor (Voigt notation) and stress &
   * strain tensors.
   * Overload for 2D Plane Strain.
   * Input the in-plane strain and the out-of-plane component.
   * Return in-plane stress and out-of-plane stress component.
   * If not providing strain_zz, default to 0
   *
   * @note This function internally converts to Voigt notation.
   */
  template <typename Tstiff, typename Tstrain, typename Tstress>
  static inline DEAL_II_ALWAYS_INLINE void
  compute_stress(const VoigtMatrix<Tstiff> &elasticity_tensor,
                 const MechTensor<Tstrain> &strain,
                 MechTensor<Tstress>       &stress,
                 Tstress                   &stress_zz)
  {
    compute_stress(elasticity_tensor, strain, Tstrain(0.0), stress, stress_zz);
  }

  /**
   * @brief Compute the stress with a elasticity tensor and stress & strain tensors.
   * Overload for 2D Plane Strain.
   * Input the in-plane strain and the out-of-plane component.
   * Return stress in Voigt notation.
   *
   * @note This function internally converts to Voigt notation.
   */
  template <typename Tstiff, typename Tstrain, typename Tstress>
  static inline DEAL_II_ALWAYS_INLINE void
  compute_stress(const VoigtMatrix<Tstiff> &elasticity_tensor,
                 const MechTensor<Tstrain> &strain,
                 const Tstrain             &strain_zz,
                 VoigtVector<Tstress>      &stress)
  {
    VoigtVector<Tstrain> epsilon = strain_to_voigt(strain, strain_zz);
    compute_stress(elasticity_tensor, epsilon, stress);
  }

  /**
   * @brief Compute the stress with a elasticity tensor and stress & strain tensors.
   * Overload for 2D Plane Strain.
   * Input the in-plane strain and the out-of-plane component.
   * Return stress in Voigt notation.
   * If not providing strain_zz, default to 0
   *
   * @note This function internally converts to Voigt notation.
   */
  template <typename Tstiff, typename Tstrain, typename Tstress>
  static inline DEAL_II_ALWAYS_INLINE void
  compute_stress(const VoigtMatrix<Tstiff> &elasticity_tensor,
                 const MechTensor<Tstrain> &strain,
                 VoigtVector<Tstress>      &stress)
  {
    compute_stress(elasticity_tensor, strain, Tstrain(0.0), stress);
  }

  // ----------------------------------------
  /**
   * @brief Isotropic stiffness matrix.
   */
  // TODO: should we use DEAL_II_ALWAYS_INLINE
  template <typename T = double>
  static inline VoigtMatrix<T>
  stiffness_isotropic(const T E, const T nu)
  {
    Assert(E > T(0.0),
           dealii::ExcMessage("Invalid isotropic elastic constants: "
                              "Young's modulus E must be positive."));
    Assert(nu > T(-1.0) && nu < T(0.5),
           dealii::ExcMessage("Invalid isotropic elastic constants: "
                              "Poisson's ratio must be in range -1 < nu < 0.5"));

    VoigtMatrix<T> stiffness;

    const T G = E / (T(2.0) * (T(1.0) + nu));

#ifdef DEBUG
    // Warning for parameter close to incompressible
    if (std::fabs(T(1.0) - T(2.0) * nu) <= tolerance<T>)
      {
        Logger::instance()
          << LogFormatter::warning(
               "WARNING: Isotropic elastic constants are nearly singular.")
          << std::endl;
      }
#endif
    // 11, 22, 33, 12
    const T lambda  = (nu * E) / ((T(1.0) + nu) * (T(1.0) - T(2.0) * nu));
    stiffness[0][0] = stiffness[1][1] = stiffness[2][2] = lambda + T(2.0) * G;
    stiffness[0][1] = stiffness[1][0] = lambda;
    stiffness[0][2] = stiffness[2][0] = lambda;
    stiffness[1][2] = stiffness[2][1] = lambda;
    stiffness[3][3]                   = G;

    return stiffness;
  }

  /**
   * @brief Orthotropic stiffness matrix: Overload for Plane Strain.
   */
  template <typename T>
  static inline VoigtMatrix<T>
  stiffness_orthotropic(const T E1,
                        const T E2,
                        const T E3,
                        const T nu12,
                        const T nu13,
                        const T nu23,
                        const T G12)
  {
    Assert(E1 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: E1 must be positive."));

    Assert(E2 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: E2 must be positive."));

    Assert(E3 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: E3 must be positive."));

    Assert(G12 > T(0.0),
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: G12 must be positive."));

    VoigtMatrix<T> stiffness;

    const T nu21 = nu12 * (E2 / E1);
    const T nu31 = nu13 * (E3 / E1);
    const T nu32 = nu23 * (E3 / E2);

    Assert(T(1.0) > nu12 * nu21,
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: 1 - nu12*nu21 must be positive."));

    Assert(T(1.0) > nu13 * nu31,
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: 1 - nu13*nu31 must be positive."));

    Assert(T(1.0) > nu23 * nu32,
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: 1 - nu23*nu32 must be positive."));

    const T delta = T(1.0) - (nu12 * nu21) - (nu23 * nu32) - (nu13 * nu31) -
                    (T(2.0) * nu12 * nu23 * nu31);

    Assert(delta > 0.0,
           dealii::ExcMessage(
             "Invalid orthotropic elastic constants: the determinant must be positive."));
#ifdef DEBUG
    // Warning for nearly singular
    if (delta <= tolerance<T>)
      {
        Logger::instance()
          << LogFormatter::warning(
               "WARNING: Orthotropic elastic constants are nearly singular.")
          << std::endl;
      }
#endif

    const T inv_delta = 1.0 / delta;

    stiffness[0][0] = E1 * (T(1.0) - nu23 * nu32) * inv_delta;
    stiffness[1][1] = E2 * (T(1.0) - nu13 * nu31) * inv_delta;
    stiffness[2][2] = E3 * (T(1.0) - nu12 * nu21) * inv_delta;
    stiffness[0][1] = stiffness[1][0] = E1 * (nu21 + nu31 * nu23) * inv_delta;
    stiffness[0][2] = stiffness[2][0] = E1 * (nu31 + nu21 * nu32) * inv_delta;
    stiffness[1][2] = stiffness[2][1] = E2 * (nu32 + nu12 * nu31) * inv_delta;
    stiffness[3][3]                   = G12;

    return stiffness;
  }

  // -------------------------------------------------
  /**
   * @brief von Mises stress (Input is in Voigt notation).
   */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE T
  stress_mises(const VoigtVector<T> &stress)
  {
    T stress_m2 = T(0);

    const T &sigma_xx = stress[0];
    const T &sigma_yy = stress[1];
    const T &sigma_zz = stress[2];
    const T &sigma_xy = stress[3];

    const T d_xy = sigma_xx - sigma_yy;
    const T d_yz = sigma_yy - sigma_zz;
    const T d_zx = sigma_zz - sigma_xx;

    stress_m2 =
      T(0.5) * (d_xy * d_xy + d_yz * d_yz + d_zx * d_zx) + T(3) * sigma_xy * sigma_xy;

    return std::sqrt(stress_m2);
  }

  //-----------------------------------------------------
  /**
 * @brief Principal stress (Input is in Voigt notation).
   @note For plane strain, returns the principal stresses of the in-plane 2x2 stress
 tensor.
 */
  template <typename T = double>
  static inline DEAL_II_ALWAYS_INLINE dealii::Tensor<1, dim, T>
                                      stress_principal(const VoigtVector<T> &stress)
  {
    dealii::Tensor<1, dim, T> stress_p {};

    constexpr unsigned int idx = 3;

    const T avg = (stress[0] + stress[1]) * T(0.5);
    const T dif = (stress[0] - stress[1]) * T(0.5);
    const T rad = std::sqrt(dif * dif + stress[idx] * stress[idx]);

    stress_p[0] = avg + rad;
    stress_p[1] = avg - rad;

    return stress_p;
  }

  // -------------------------------------------------------

  /**
   * @brief Extract 4x4 plane strain stiffness from 6x6 3D stiffness.
   * @note This restricted form assumes the xz and yz elastic shear components are zero,
   * which may not work for general anisotropy and eigenstrain.
   */
  template <typename T>
  static inline DEAL_II_ALWAYS_INLINE VoigtMatrix<T>
  extract_plane_strain_stiffness(const ThreeDimensionalStiffness<T> &stiffness_3d)
  {
    VoigtMatrix<T> stiffness;

    constexpr std::array<unsigned int, 4> index = {0, 1, 2, 5};

    for (unsigned int i = 0; i < 4; ++i)
      {
        for (unsigned int j = 0; j < 4; ++j)
          {
            stiffness[i][j] = stiffness_3d[index[i]][index[j]];
          }
      }

    return stiffness;
  }
};

template <typename T = double>
using PlaneStrainVector = PlaneStrain::VoigtVector<T>;

template <typename T = double>
using PlaneStrainStiffness = PlaneStrain::VoigtMatrix<T>;

namespace OldMechanics
{
  /**
   * -------------------------------------------------------------
   * The following are provided just for backward compatibility.
   * It is preferred to use the functions above when calling
   * compute_stress.
   * -------------------------------------------------------------
   */

  /**
   * @brief Voigt notation index range.
   * This is used for backward compatibility.
   */
  template <unsigned int dim>
  [[deprecated("Use voigt_size<dim, StressState> instead.")]] constexpr unsigned int
    voigt_tensor_size = (2 * dim) - 1 + (dim / 3);

  /**
   * @brief Compute the stress with a given displacement and elasticity tensor. This
   * assumes that the provided parameters are in Voigt notation.
   */
  template <unsigned int dim, typename T>
  static inline DEAL_II_ALWAYS_INLINE void
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
  static inline DEAL_II_ALWAYS_INLINE void
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
}; // namespace OldMechanics

PRISMS_PF_END_NAMESPACE
// NOLINTEND(readability-identifier-naming,readability-magic-numbers,cppcoreguidelines-avoid-magic-numbers,readability-identifier-length)
