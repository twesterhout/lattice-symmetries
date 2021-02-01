cmake_minimum_required(VERSION 3.15)

project(${CMAKE_PROJECT_NAME}Tests LANGUAGES C CXX)

set(test_sources 
    src/test_main.cpp
    src/test_api.cpp
    # src/test_permutation.cpp
    # src/test_network.cpp
    # src/test_cache.cpp
)
add_executable(${PROJECT_NAME} ${test_sources})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

if (${CMAKE_PROJECT_NAME}_ENABLE_CODE_COVERAGE)
  target_compile_options(${CMAKE_PROJECT_NAME} PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
  target_link_options(${CMAKE_PROJECT_NAME} PUBLIC -fprofile-arcs -ftest-coverage)
endif()

find_package(Catch2 QUIET)
if(NOT Catch2_FOUND)
  message(STATUS "Catch2 not found. Downloading")
  include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.0
  ) 

  FetchContent_GetProperties(Catch2)
  if(NOT catch2_POPULATED)
    FetchContent_Populate(Catch2)
    add_subdirectory(${catch2_SOURCE_DIR} ${catch2_BINARY_DIR})
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/contrib")
  endif()

  FetchContent_MakeAvailable(Catch2)
  message(STATUS "Catch2 not found. Downloading - done")
endif()

target_link_libraries(
  ${PROJECT_NAME} 
  PUBLIC
    Catch2::Catch2
    lattice_symmetries
  PRIVATE
    lattice_symmetries_search_sorted
)
# Header files in src/ are not part of the external interface, but we still
# need to test them.
target_include_directories(
    ${PROJECT_NAME}
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/third_party
)

include(Catch)
catch_discover_tests(${PROJECT_NAME})
