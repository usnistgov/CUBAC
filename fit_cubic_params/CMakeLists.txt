cmake_minimum_required(VERSION 2.8)

project(example)

# Add CoolProp dependency
set(COOLPROP_STATIC_LIBRARY true)
add_subdirectory("${CMAKE_SOURCE_DIR}/../externals/CoolProp" "CoolProp")

# Add simple example
add_executable(Main "${CMAKE_SOURCE_DIR}/cubac.cpp")
target_link_libraries (Main CoolProp)

# Add pybind11 wrapper
add_subdirectory("${CMAKE_SOURCE_DIR}/../externals/pybind11" "pybind11")
pybind11_add_module(PureFluid "${CMAKE_SOURCE_DIR}/cubac.cpp")
target_link_libraries(PureFluid CoolProp)
target_compile_definitions(PureFluid PUBLIC -DPYBIND11)
if ( MSVC )
  set_target_properties( PureFluid PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR})
endif ()