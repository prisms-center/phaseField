#
# Find libassert and run some checks
#

# With libassert we always want to vendor it

# Add the external projects
prisms_pf_add_external_project(
  libassert
  https://github.com/jeremy-rifkin/libassert.git
  v2.2.1
)

# Create the libraries
prisms_pf_add_external_library(libassert libassert)

# Create install rules
prisms_pf_install_external_library(libassert)

# Add libassert to the Release and Debug lists
prisms_pf_add_dependency_target(imported_libassert_debug DEBUG PUBLIC)
prisms_pf_add_dependency_target(imported_libassert RELEASE PUBLIC)
