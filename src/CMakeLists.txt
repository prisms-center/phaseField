# This file collects the source files for PRISMS-PF
set(PRISMS_PF_SOURCE_FILES)
foreach(TARGET IN ITEMS
  core
  field_input
  grains
  nucleation
  utilities
)
  add_subdirectory(${TARGET})
endforeach()
set(PRISMS_PF_SOURCE_FILES ${PRISMS_PF_SOURCE_FILES} PARENT_SCOPE)