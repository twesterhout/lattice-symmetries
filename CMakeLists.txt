cmake_minimum_required(VERSION 3.15)

project(
  LatticeSymmetries
  VERSION 0.1.1
  LANGUAGES C CXX
)

#
# Prevent building in the source directory
#
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.")
endif()

#
# Project settings
#
option(${PROJECT_NAME}_WARNINGS_AS_ERRORS "Treat compiler warnings as errors." OFF)
option(${PROJECT_NAME}_ENABLE_UNIT_TESTING "Enable unit tests for the project." ON)
option(${PROJECT_NAME}_ENABLE_CLANG_TIDY "Enable static analysis with clang-tidy." OFF)
option(${PROJECT_NAME}_ENABLE_CPPCHECK "Enable static analysis with Cppcheck." OFF)
option(${PROJECT_NAME}_ENABLE_CODE_COVERAGE "Enable code coverage through GCC." OFF)
option(${PROJECT_NAME}_LINK_STDLIB_STATICALLY "Enable static linking of libstdc++ and libgcc on Linux." OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS OFF)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

#
# Set the project standard and warnings
#
include(cmake/CompilerWarnings.cmake)

#
# Sanitizers
#
include(cmake/StaticAnalyzers.cmake)

#
# Library sources
#
set(LatticeSymmetries_sources
    src/permutation.cpp
    src/network.cpp
    src/symmetry.cpp
    src/group.cpp
    src/cache.cpp
    src/basis.cpp
    src/operator.cpp
    src/error_handling.cpp)

set(LatticeSymmetries_headers
    include/lattice_symmetries/lattice_symmetries.h
    src/permutation.hpp
    src/network.hpp
    src/symmetry.hpp
    src/cache.hpp
    src/operator.hpp
    src/error_handling.hpp)

set(LatticeSymmetries_all_files
    ${LatticeSymmetries_sources}
    ${LatticeSymmetries_headers}
    src/kernels.hpp
    src/kernels.cpp)


#
# Add a target for formating the project using `clang-format` (i.e: cmake --build build --target clang-format)
#
function(add_clang_format_target)
    if(NOT ${PROJECT_NAME}_CLANG_FORMAT_BINARY)
        find_program(${PROJECT_NAME}_CLANG_FORMAT_BINARY
            NAMES clang-format
                  clang-format-11
                  clang-format-10)
    endif()
    if(${PROJECT_NAME}_CLANG_FORMAT_BINARY)
        add_custom_target(clang-format
            COMMAND ${${PROJECT_NAME}_CLANG_FORMAT_BINARY}
            -i ${CMAKE_CURRENT_LIST_DIR}/${LatticeSymmetries_all_files})
    else()
        message(WARNING "`clang-format` target is disabled because clang-format binary was not found")
    endif()
endfunction()

add_clang_format_target()

function(disable_rtti_and_exceptions _target)
    target_compile_options(${_target} PUBLIC -fno-rtti -fno-exceptions)
endfunction()

#
# CPU kernels for various architectures
#
foreach(_target sse2 avx avx2)
    add_library(lattice_symmetries_kernels_${_target}
        OBJECT src/kernels.hpp
               src/kernels.cpp)
    target_include_directories(lattice_symmetries_kernels_${_target}
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    )
    target_include_directories(lattice_symmetries_kernels_${_target}
      SYSTEM PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party
    )
    target_compile_features(lattice_symmetries_kernels_${_target} PRIVATE cxx_std_17)
    set_project_warnings(lattice_symmetries_kernels_${_target})
    disable_rtti_and_exceptions(lattice_symmetries_kernels_${_target})
    set_property(TARGET lattice_symmetries_kernels_${_target} PROPERTY
        POSITION_INDEPENDENT_CODE ON)
endforeach()
target_compile_definitions(lattice_symmetries_kernels_sse2 PRIVATE LATTICE_SYMMETRIES_ADD_DISPATCH_CODE=1)
target_compile_options(lattice_symmetries_kernels_sse2 PRIVATE -m64 -msse2)
target_compile_options(lattice_symmetries_kernels_avx PRIVATE -m64 -mavx)
target_compile_options(lattice_symmetries_kernels_avx2 PRIVATE -m64 -mavx2 -mfma)

#
# Create library
#
add_library(lattice_symmetries
    ${LatticeSymmetries_headers}
    ${LatticeSymmetries_sources}
)
target_compile_features(lattice_symmetries PRIVATE cxx_std_17)
set_project_warnings(lattice_symmetries)
disable_rtti_and_exceptions(lattice_symmetries)
set_property(TARGET lattice_symmetries PROPERTY POSITION_INDEPENDENT_CODE ON)
target_sources(lattice_symmetries
    PRIVATE
        $<TARGET_OBJECTS:lattice_symmetries_kernels_sse2>
        $<TARGET_OBJECTS:lattice_symmetries_kernels_avx>
        $<TARGET_OBJECTS:lattice_symmetries_kernels_avx2>)
if(APPLE)
    message(WARNING "C++ standard library will be linked dynamically")
else()
    if(${PROJECT_NAME}_LINK_STDLIB_STATICALLY)
        message(STATUS "Linking libstdc++ and libgcc statically")
        target_link_options(lattice_symmetries PRIVATE -static-libstdc++ -static-libgcc)
    endif()
endif()

# Allow usage of header files in the `src` directory, but only for utilities
target_include_directories(
  lattice_symmetries
  PUBLIC 
    $<INSTALL_INTERFACE:include>    
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_include_directories(
  lattice_symmetries
  SYSTEM PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party
)

#
# Dependencies 
#
find_package(OpenMP REQUIRED)

target_link_libraries(
  lattice_symmetries
  PUBLIC 
    OpenMP::OpenMP_CXX
)


#
# Provide alias to library for 
#
add_library(LatticeSymmetries::LatticeSymmetries ALIAS lattice_symmetries)

#
# Install library for easy downstream inclusion
#
include(GNUInstallDirs)

set(LATTICE_SYMMETRIES_PC_PREFIX ${CMAKE_INSTALL_PREFIX})
set(LATTICE_SYMMETRIES_PC_VERSION ${PROJECT_VERSION})
set(LATTICE_SYMMETRIES_PC_LIBS "-L\${libdir} -Wl,-rpath=\${libdir} -llattice_symmetries")
set(LATTICE_SYMMETRIES_PC_LIBS_PRIVATE "")
configure_file(
    "${CMAKE_CURRENT_LIST_DIR}/cmake/lattice_symmetries.pc.in"
    "${CMAKE_BINARY_DIR}/lattice_symmetries.pc"
    @ONLY
)

install(
  TARGETS
    lattice_symmetries
)
install(
  FILES
    "${CMAKE_BINARY_DIR}/lattice_symmetries.pc"
  DESTINATION
    lib/pkgconfig
)
install(
  DIRECTORY
    include/
  DESTINATION
    include
)

#
# Unit testing
#
if(${PROJECT_NAME}_ENABLE_UNIT_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

#
# Benchmarks
#
# add_subdirectory(benchmark)
