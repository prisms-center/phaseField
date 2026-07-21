#
# Find Catch2 and run some checks
#

# With Catch2 we always want to vendor it
# NOTE: We use FetchContent here for Catch2 to make
# things easier. With the unit tests we don't have
# to worry about DebugRelease anyway.

include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.13.0
  SYSTEM
)
FetchContent_MakeAvailable(Catch2)
