cmake_minimum_required(VERSION 3.30.0)
project(QKD_Winnow_V VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM.cmake")
if(NOT(EXISTS ${CPM_DOWNLOAD_LOCATION}))
    message(STATUS "Downloading CPM.cmake")
    file(DOWNLOAD https://github.com/cpm-cmake/CPM.cmake/releases/latest/download/CPM.cmake ${CPM_DOWNLOAD_LOCATION})
endif()
include(${CPM_DOWNLOAD_LOCATION})

CPMAddPackage(
    NAME BS_thread_pool
    GITHUB_REPOSITORY bshoshany/thread-pool
    VERSION 4.1.0)
add_library(BS_thread_pool INTERFACE)
target_include_directories(BS_thread_pool INTERFACE ${BS_thread_pool_SOURCE_DIR}/include)

CPMAddPackage(
    NAME fmt
    GIT_TAG 11.0.2
    GITHUB_REPOSITORY fmtlib/fmt)

CPMAddPackage(
    NAME indicators
    GITHUB_REPOSITORY p-ranav/indicators
    VERSION 2.3)

CPMAddPackage(
    NAME xoshiro_cpp
    GITHUB_REPOSITORY Reputeless/Xoshiro-cpp
    VERSION 1.1)
include_directories(${xoshiro_cpp_SOURCE_DIR})

CPMAddPackage(
    NAME nlohmann_json
    VERSION 3.11.3
    GITHUB_REPOSITORY nlohmann/json
    OPTIONS
        "JSON_BuildTests OFF")

option(USE_CURRENT_DIR "Use current directory where executable file is located instead of SOURCE_DIR to read and write files" OFF)
if(USE_CURRENT_DIR)
    add_definitions(-DUSE_CURRENT_DIR)
endif()
add_definitions(-DSOURCE_DIR="${CMAKE_SOURCE_DIR}")

add_executable(QKD_Winnow_V 
    src/main.cpp
    src/config.hpp
    src/config.cpp
    src/utils.hpp
    src/utils.cpp
    src/simulation.hpp
    src/simulation.cpp
    src/bit_array_operations.hpp
    src/bit_array_operations.cpp
    src/winnow_algorithm.hpp
    src/winnow_algorithm.cpp
    )
target_include_directories(QKD_Winnow_V PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(QKD_Winnow_V PUBLIC BS_thread_pool fmt indicators nlohmann_json)

option(BUILD_EXAMPLES "Build examples" ON)

if(BUILD_EXAMPLES)
    add_executable(QKD_Winnow_V_Example 
    example/winnow_example.cpp
    src/config.hpp
    src/config.cpp
    src/utils.hpp
    src/utils.cpp
    src/bit_array_operations.hpp
    src/bit_array_operations.cpp
    src/winnow_algorithm.hpp
    src/winnow_algorithm.cpp
    )
    target_include_directories(QKD_Winnow_V_Example PRIVATE ${CMAKE_SOURCE_DIR}/src)
    target_link_libraries(QKD_Winnow_V_Example fmt nlohmann_json)
endif()