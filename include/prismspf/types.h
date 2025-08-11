#ifndef types_h
#define types_h

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

namespace types
{
  /**
   * \brief Type for field indices.
   */
  using index = unsigned int;

} // namespace types

namespace numbers
{
  /**
   * \brief Invalid field index.
   */
  static const prisms::types::index invalid_index = static_cast<prisms::types::index>(-1);

} // namespace numbers

PRISMS_PF_END_NAMESPACE

#endif