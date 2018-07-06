cmake_minimum_required(VERSION 3.10)
project(tba-solver)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

# set build type to Release
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

include_directories(${CMAKE_SOURCE_DIR}/src)

# examples
add_custom_target(examples)
file(GLOB ex_src ${CMAKE_SOURCE_DIR}/examples/*/*.cpp)
foreach(ex ${ex_src})
  get_filename_component(tag ${ex} NAME_WE)
  add_executable(${tag} ${ex})
  install(TARGETS ${tag} DESTINATION ${CMAKE_SOURCE_DIR}/examples)
  add_dependencies(examples ${tag})
endforeach()
