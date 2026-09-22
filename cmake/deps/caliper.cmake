#
# Find Caliper and run some checks
#

# TODO: Add config variables that tell us about caliper here

find_package(
  caliper
  REQUIRED
  HINTS
    ${CALIPER_DIR}
    $ENV{CALIPER_DIR}
)
