#
# Find Boost.JSON
#

if(NOT Boost_ROOT AND DEFINED ENV{BOOST_ROOT})
  set(Boost_ROOT "$ENV{BOOST_ROOT}")
endif()

if(Boost_ROOT)
  file(GLOB _boost_config_files "${Boost_ROOT}/lib/cmake/Boost-*/BoostConfig.cmake")
  if(_boost_config_files)
    list(GET _boost_config_files -1 _boost_config_file)
    get_filename_component(Boost_DIR "${_boost_config_file}" DIRECTORY)
  endif()
endif()

find_package(Boost CONFIG REQUIRED COMPONENTS json)

prisms_pf_add_dependency_targets(DEBUG PUBLIC Boost::json)
prisms_pf_add_dependency_targets(RELEASE PUBLIC Boost::json)
