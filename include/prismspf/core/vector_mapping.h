#pragma once

#include <deal.II/lac/block_vector.h>

#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <cstdint>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Information about the field.
 */
struct FieldInfo
{
  /**
   * @brief Tensor rank of the field.
   *
   * Currently, only scalar and vectors are supported.
   */
  enum class TensorRank : std::uint8_t
  {
    Undefined,
    Scalar,
    Vector
  };

  /**
   * @brief Type of PDE that is being solved.
   */
  enum class PDEType : std::uint8_t
  {
    Undefined,
    Constant,
    ExplicitTimeDependent,
    ImplicitTimeDependent,
    TimeIndependent,
    Auxiliary
  };

  /**
   * @brief Internal classification of PDE types.
   *
   * There are several different types of solves that are possible.
   *
   * For PDEType::ExplicitTimeDependent, all fields of that type can be solved
   * concurrently. This also applies for postprocessed fields of that type.
   *
   * For PDEType::Constant, those fields are never updated.
   *
   * For the rest, we have several classifications that can either be solved sequentially
   * of concurrently.
   */
  enum class SolveType : std::uint8_t
  {
    Undefined,
    ConcurrentConstant,
    ConcurrentExplicit,
    ConcurrentAuxiliary,
    SequentialAuxiliary,
    SequentialLinear,
    SequentialSelfNonlinear,
    SequentialCoNonlinear
  };

  /**
   * @brief Internal classification of dependency type flags.
   */
  enum class DependencyTypeFlags : std::uint8_t
  {
    None = 0x00,

    New     = 0x01,
    Current = 0x02,

    Change = 0x04,

    OldOne   = 0x08,
    OldTwo   = 0x10,
    OldThree = 0x20,
    OldFour  = 0x40
  };

  friend DependencyTypeFlags
  operator|(const DependencyTypeFlags flag_1, const DependencyTypeFlags flag_2)
  {
    return static_cast<DependencyTypeFlags>(static_cast<std::uint8_t>(flag_1) |
                                            static_cast<std::uint8_t>(flag_2));
  }

  friend DependencyTypeFlags &
  operator|=(DependencyTypeFlags &flag_1, const DependencyTypeFlags flag_2)
  {
    flag_1 = flag_1 | flag_2;
    return flag_1;
  }

  friend DependencyTypeFlags
  operator&(const DependencyTypeFlags flag_1, const DependencyTypeFlags flag_2)
  {
    return static_cast<DependencyTypeFlags>(static_cast<std::uint8_t>(flag_1) &
                                            static_cast<std::uint8_t>(flag_2));
  }

  friend DependencyTypeFlags &
  operator&=(DependencyTypeFlags &flag_1, const DependencyTypeFlags flag_2)
  {
    flag_1 = flag_1 & flag_2;
    return flag_1;
  }

  bool                is_postprocess        = false;
  bool                is_nucleation_rate    = false;
  TensorRank          tensor_rank           = TensorRank::Undefined;
  PDEType             pde_type              = PDEType::Undefined;
  SolveType           solve_type            = SolveType::Undefined;
  DependencyTypeFlags dependency_type_flags = DependencyTypeFlags::None;
  unsigned int        global_index          = Numbers::invalid_index;
  unsigned int        solve_block_index     = Numbers::invalid_index;
};

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
