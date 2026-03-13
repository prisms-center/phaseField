#
# Decompose version into major, minor, and subminor
#

function(version_decomposition VERSION)
  # Replace any `.` with a `;` so we get a list
  string(REPLACE "." ";" VERSION_LIST ${VERSION})

  list(GET VERSION_LIST 0 MAJOR)
  list(GET VERSION_LIST 1 MINOR)
  list(GET VERSION_LIST 2 SUBMINOR)

  set(PRISMS_PF_VERSION_MAJOR ${MAJOR} PARENT_SCOPE)
  set(PRISMS_PF_VERSION_MINOR ${MINOR} PARENT_SCOPE)
  set(PRISMS_PF_VERSION_SUBMINOR ${SUBMINOR} PARENT_SCOPE)
endfunction()
