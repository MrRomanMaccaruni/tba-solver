cmake_minimum_required(VERSION 3.10)
project(tba-solver)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

include_directories(${CMAKE_SOURCE_DIR}/src)

# examples
add_custom_target(examples)
file(GLOB ex_src ${CMAKE_SOURCE_DIR}/examples/*lee*/*lee*.cpp)
foreach(ex ${ex_src})
  get_filename_component(tag ${ex} NAME_WE)
  add_executable(${tag} ${ex})
  install(TARGETS ${tag} DESTINATION ${CMAKE_SOURCE_DIR}/examples)
  add_dependencies(examples ${tag})
endforeach()
