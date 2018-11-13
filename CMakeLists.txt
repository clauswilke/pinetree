cmake_minimum_required(VERSION 3.1.0)

project(pinetree)

# Modern C++ woooo!
set(CMAKE_CXX_STANDARD 11)

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
install(TARGETS core DESTINATION src/pinetree)
