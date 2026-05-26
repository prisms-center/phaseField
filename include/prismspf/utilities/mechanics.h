#pragma once

#include <deal.II/base/tensor.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

namespace Mechanics
{

  /**
   * @brief Voigt notation index range
   */
  template <unsigned int dim>
  constexpr unsigned int voigt_tensor_size = (2 * dim) - 1 + (dim / 3);

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
