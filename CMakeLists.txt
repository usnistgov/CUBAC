cmake_minimum_required (VERSION 2.8)
project (CUBAC)

# All source files to be compiled shall have the .cpp file extension and reside in the src directory
file(GLOB_RECURSE APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# All global includes shall reside in the include directory
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/rapidjson/include")

add_executable(Poling ${CMAKE_CURRENT_SOURCE_DIR}/examples/Poling.cpp ${APP_SOURCES})
add_executable(methane_water ${CMAKE_CURRENT_SOURCE_DIR}/examples/methane_water.cpp ${APP_SOURCES})
add_executable(speed_test ${CMAKE_CURRENT_SOURCE_DIR}/examples/speed_test.cpp ${APP_SOURCES})
add_executable(VTPR ${CMAKE_CURRENT_SOURCE_DIR}/examples/VTPR.cpp ${APP_SOURCES})