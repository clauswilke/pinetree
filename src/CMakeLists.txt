cmake_minimum_required (VERSION 2.6)
project (Pysinthe)
set (CMAKE_CXX_STANDARD 11)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake.modules/")
enable_testing()
                                                                                                 
file(GLOB SOURCE_FILES "*.cpp" "./tests/*.cpp")                                                          
add_executable(${PROJECT_NAME} ${SOURCE_FILES})                                          
include(ParseAndAddCatchTests)                                                           
ParseAndAddCatchTests(${PROJECT_NAME})  
