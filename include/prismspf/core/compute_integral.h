// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef compute_integral_h
#define compute_integral_h

#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/config.h>
#include <prismspf/core/element_volume.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Compute the integral of a given field.
 */
template <int dim, int degree, typename number>
class computeIntegral
{
public:
  /**
   * \brief Constructor.
   */
  explicit computeIntegral(const elementVolume<dim, degree, number> &_element_volume);

  /**
   * \brief Destructor.
   */
  ~computeIntegral() = default;

  /**
   * \brief Initialize.
   */
  void
  initialize(
    std::shared_ptr<dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>
      _data);

  /**
   * \brief Compute the integral
   */
  void
  compute_integral();

private:
  /**
   * \brief Element volumes
   */
  const elementVolume<dim, degree, double> &element_volume;

  /**
   * \brief Matrix-free object.
   */
  std::shared_ptr<dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>> data;
};

template <int dim, int degree, typename number>
computeIntegral<dim, degree, number>::computeIntegral(
  const elementVolume<dim, degree, number> &_element_volume)
  : element_volume(_element_volume)
{}

template <int dim, int degree, typename number>
void
computeIntegral<dim, degree, number>::initialize(
  std::shared_ptr<dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>> _data)
{
  data = _data;
}

template <int dim, int degree, typename number>
void
computeIntegral<dim, degree, number>::compute_integral()
{}

PRISMS_PF_END_NAMESPACE

#endif