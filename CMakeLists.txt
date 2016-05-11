cmake_minimum_required (VERSION 2.8)
project (CUBAC)

# All source files to be compiled shall have the .cpp file extension and reside in the src directory
file(GLOB_RECURSE APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# All global includes shall reside in the include directory
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/rapidjson/include")

add_executable(Poling ${CMAKE_CURRENT_SOURCE_DIR}/examples/Poling.cpp ${APP_SOURCES})
