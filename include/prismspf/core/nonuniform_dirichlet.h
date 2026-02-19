// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

template <unsigned int dim, unsigned int degree, typename number>
class PDEOperatorBase;

/**
 * @brief Function for user-implemented nonuniform dirichlet boundary condition.
 */
template <unsigned int dim, unsigned int degree, typename number>
class NonuniformDirichlet : public dealii::Function<dim, number>
{
public:
  /**
   * @brief Constructor.
   */
  NonuniformDirichlet(unsigned int                                _index,
                      unsigned int                                _boundary_id,
                      const PDEOperatorBase<dim, degree, number> &_pde_operator,
                      unsigned int                                spacedim);

  // NOLINTBEGIN(readability-identifier-length, readability-avoid-const-params-in-decls)

  /**
   * @brief Scalar value.
   */
  number
  value(const dealii::Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * @brief Vector value.
   */
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<number> &value) const override;

  // NOLINTEND(readability-identifier-length, readability-avoid-const-params-in-decls)

private:
  unsigned int index;

  unsigned int boundary_id;

  std::shared_ptr<const PDEOperatorBase<dim, degree, number>> pde_operator;
};

PRISMS_PF_END_NAMESPACE
