#pragma once

#include <deal.II/lac/block_vector.h>

#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * Due the generalized nature of PRISMS-PF and how we handle the solving of PDEs on the
 * backend, we end up with a slew of vectors that have to be managed and updated as we
 * solve fields. While we can stored these in hierarchical data structure to perserve the
 * logic and intuition this becomes unwieldy on the backend and is less performant than
 * storing them in a flat structure or vector blocks.
 *
 * This class has a unified description of these mappings so that everywhere else in the
 * code we can simply use 1D vector indices.
 */
template <typename Number>
class VectorMappings
{
public:
  /**
   * @brief Constructor.
   */
  explicit VectorMappings(const std::vector<VariableAttributes> &variable_attributes);

private:
  /**
   * @brief Number of blocks.
   *
   * The number of blocks is determined by the dependencies of each field upon
   * construction of this class.
   */
  unsigned int n_blocks = 0;

  /**
   * @brief Block vector of all neccessary vectors for the PDE problem.
   */
  dealii::BlockVector<Number> block_vector;

  /**
   *
   */
};

PRISMS_PF_END_NAMESPACE
