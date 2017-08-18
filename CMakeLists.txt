cmake_minimum_required(VERSION 2.8.12)

project(pinetree)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/lib/cmake.modules/")

# For VSCode debugger support
include(CMakeToolsHelpers OPTIONAL)
# Required for CTest to function properly
include(Dart)
# Automatically parse Catch unit tests
include(ParseAndAddCatchTests)
# Set up for Catch unit testing
enable_testing()

# Modern C++ woooo!
set(CMAKE_CXX_STANDARD 11)

# Tell Cmake that headers are also in source directory
include_directories(${PROJECT_NAME})
include_directories(tests)
include_directories(lib/catch/include)

# Pull in source and test files
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/${PROJECT_NAME}")
set(TEST_DIR "${CMAKE_SOURCE_DIR}/tests")
set(SOURCES 
    "${SOURCE_DIR}/feature.cpp"
    "${SOURCE_DIR}/polymer.cpp"
    "${SOURCE_DIR}/choices.cpp"
    "${SOURCE_DIR}/tracker.cpp"
    "${SOURCE_DIR}/simulation.cpp")
set(TESTS ${SOURCES}
    "${TEST_DIR}/test_main.cpp"
    "${TEST_DIR}/test_feature.cpp"
    "${TEST_DIR}/test_polymer.cpp"
    "${TEST_DIR}/test_choices.cpp"
    "${TEST_DIR}/test_signals.cpp"
    "${TEST_DIR}/test_simulation.cpp")

# Generate a test executable
add_executable("${PROJECT_NAME}_test" ${TESTS})

# Automatically parse Catch unit tests
ParseAndAddCatchTests("${PROJECT_NAME}_test")

# Python bindings!
add_subdirectory(${CMAKE_SOURCE_DIR}/lib/pybind11)
pybind11_add_module(${PROJECT_NAME} "${SOURCE_DIR}/python_bindings.cpp" ${SOURCES})
