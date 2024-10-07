cmake_minimum_required(VERSION 3.9.0)

project(pinetree)

# Modern C++ woooo!
set(CMAKE_CXX_STANDARD 11)

# uncomment for OpenMP suport

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0 -g")

# Tell Cmake that headers are also in source directory
include_directories(src/${PROJECT_NAME})
include_directories(tests)
include_directories(lib)

# Pull in source and test files
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/${PROJECT_NAME}")
set(SOURCES 
    "${SOURCE_DIR}/feature.cpp"
    "${SOURCE_DIR}/polymer.cpp"
    "${SOURCE_DIR}/choices.cpp"
    "${SOURCE_DIR}/tracker.cpp"
    "${SOURCE_DIR}/model.cpp"
    "${SOURCE_DIR}/gillespie.cpp"
    "${SOURCE_DIR}/reaction.cpp")

# Generate python module
add_subdirectory(lib/pybind11)
pybind11_add_module(core ${SOURCES} "${SOURCE_DIR}/python_bindings.cpp")
install(TARGETS core DESTINATION src/${PROJECT_NAME})

SET(TEST_DIR "tests")
SET(TESTS ${SOURCES}
    "${TEST_DIR}/test_main.cpp"
    "${TEST_DIR}/unit_tests.cpp")

# Generate a test executable
#include_directories(lib/catch/include)
add_executable("${PROJECT_NAME}_test" ${TESTS})
