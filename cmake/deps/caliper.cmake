#
# Find Caliper and run some checks
#

find_package(
  caliper
  REQUIRED
  HINTS
    ${CALIPER_DIR}
    $ENV{CALIPER_DIR}
)
