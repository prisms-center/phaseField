// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/solver_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolveContext;

/**
 * @brief This class handles the solves for fields that are constant in time (does
 * nothing that the base solver class doesn't already do).
 */
template <unsigned int dim, unsigned int degree, typename number>
using ConstantSolver = SolverBase<dim, degree, number>;

PRISMS_PF_END_NAMESPACE
