cmake_minimum_required(VERSION 3.15)

# Read version from __init__.py file
file(STRINGS "${CMAKE_CURRENT_LIST_DIR}/python/lattice_symmetries/__init__.py" _init_py_contents)
set(_version_regex "^__version__[ \t\r\n]*=[ \t\r\n]*['\"]([^'\"]*)['\"]")
# Pick the correct line
list(FILTER _init_py_contents INCLUDE REGEX "${_version_regex}")
# Extract the version
string(REGEX REPLACE ${_version_regex} "\\1" _version "${_init_py_contents}")

project(
  LatticeSymmetries
  VERSION ${_version}
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
option(${PROJECT_NAME}_ENABLE_PROFILING "Enable profiling" OFF)
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
    src/basis.cpp
    src/cache.cpp
    src/error_handling.cpp
    src/group.cpp
    src/network.cpp
    src/operator.cpp
    src/permutation.cpp
    src/symmetry.cpp
    src/cpu/binary_search.cpp
)
set(LatticeSymmetries_portable_sources
    src/error_handling.cpp
    src/group.cpp
    src/permutation.cpp
    src/symmetry.cpp
)

set(LatticeSymmetries_headers
    include/lattice_symmetries/lattice_symmetries.h
    src/basis.hpp
    src/bits.hpp
    src/cache.hpp
    src/intrusive_ptr.hpp
    src/network.hpp
    src/operator.hpp
    src/permutation.hpp
    src/symmetry.hpp
)

set(LatticeSymmetries_all_files
    ${LatticeSymmetries_sources}
    ${LatticeSymmetries_headers}
)


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
    target_compile_options(${_target} PRIVATE -fno-rtti -fno-exceptions)
endfunction()

#
# CPU kernels for various architectures
#
foreach(arch sse2 sse4 avx avx2)
    set(_local_target "lattice_symmetries_kernels_${arch}")
    add_library(${_local_target} OBJECT
        src/cpu/search_sorted.cpp
        src/cpu/benes_forward_64.cpp
        src/cpu/benes_forward_512.cpp
        src/cpu/state_info.cpp
    )
    target_include_directories(${_local_target}
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_include_directories(${_local_target}
      SYSTEM PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/version2
    )
    target_compile_features(${_local_target} PRIVATE cxx_std_17)
    target_compile_definitions(${_local_target} PRIVATE VCL_NAMESPACE=vcl_${arch})
    target_compile_definitions(${_local_target} PRIVATE OUTCOME_DISABLE_EXECINFO)
    set_project_warnings(${_local_target})
    disable_rtti_and_exceptions(${_local_target})
    set_property(TARGET ${_local_target} PROPERTY POSITION_INDEPENDENT_CODE ON)
    # target_link_libraries(${_local_target} PUBLIC OpenMP::OpenMP_CXX)
endforeach()
target_compile_definitions(lattice_symmetries_kernels_sse2 PRIVATE LATTICE_SYMMETRIES_ADD_DISPATCH_CODE=1)
target_compile_options(lattice_symmetries_kernels_sse2 PRIVATE -m64 -march=nocona)
target_compile_options(lattice_symmetries_kernels_sse4 PRIVATE -msse4.1 -msse4.2 -mpopcnt)
target_compile_options(lattice_symmetries_kernels_avx PRIVATE -mavx)
target_compile_options(lattice_symmetries_kernels_avx2 PRIVATE -mavx2 -mbmi2 -mfma)


# foreach(_target sse2 avx avx2)
#     add_library(lattice_symmetries_kernels_${_target}
#         OBJECT src/kernels.hpp
#                src/kernels.cpp
#                # src/state_info.hpp
#                # src/state_info.cpp
#     )
#     target_include_directories(lattice_symmetries_kernels_${_target}
#       PRIVATE
#         ${CMAKE_CURRENT_SOURCE_DIR}/include
#         ${CMAKE_CURRENT_SOURCE_DIR}/src
#     )
#     target_include_directories(lattice_symmetries_kernels_${_target}
#       SYSTEM PRIVATE
#         ${CMAKE_CURRENT_SOURCE_DIR}/third_party
#         ${CMAKE_CURRENT_SOURCE_DIR}/third_party/version2
#     )
#     target_compile_features(lattice_symmetries_kernels_${_target} PRIVATE cxx_std_17)
#     target_compile_definitions(lattice_symmetries_kernels_${_target} PRIVATE VCL_NAMESPACE=vcl_${_target})
#     set_project_warnings(lattice_symmetries_kernels_${_target})
#     disable_rtti_and_exceptions(lattice_symmetries_kernels_${_target})
#     set_property(TARGET lattice_symmetries_kernels_${_target} PROPERTY
#         POSITION_INDEPENDENT_CODE ON)
# endforeach()
# target_compile_definitions(lattice_symmetries_kernels_sse2 PRIVATE LATTICE_SYMMETRIES_ADD_DISPATCH_CODE=1)
# target_compile_options(lattice_symmetries_kernels_sse2 PRIVATE -m64 -msse2)
# target_compile_options(lattice_symmetries_kernels_avx PRIVATE -m64 -mavx)
# target_compile_options(lattice_symmetries_kernels_avx2 PRIVATE -m64 -mavx2 -mfma)


#
# Pure C code
#
add_library(lattice_symmetries_c_sources OBJECT
    src/error_handling.c
    src/dispatcher.c
    src/batched.c
)
set_property(TARGET lattice_symmetries_c_sources PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(lattice_symmetries_c_sources
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

#
# Create library
#
add_library(lattice_symmetries
    ${LatticeSymmetries_headers}
    ${LatticeSymmetries_sources}
)
target_link_libraries(lattice_symmetries PRIVATE lattice_symmetries_c_sources)
foreach(arch sse2 sse4 avx avx2)
    target_link_libraries(lattice_symmetries PRIVATE lattice_symmetries_kernels_${arch})
endforeach()


target_compile_features(lattice_symmetries PRIVATE cxx_std_17)
set_project_warnings(lattice_symmetries)
disable_rtti_and_exceptions(lattice_symmetries)
target_compile_definitions(lattice_symmetries PRIVATE OUTCOME_DISABLE_EXECINFO)
set_property(TARGET lattice_symmetries PROPERTY POSITION_INDEPENDENT_CODE ON)
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
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/version2
)
target_compile_definitions(lattice_symmetries PRIVATE VCL_NAMESPACE=vcl)

#
# Dependencies 
#
find_package(OpenMP REQUIRED)

target_link_libraries(
  lattice_symmetries
  PUBLIC 
    OpenMP::OpenMP_CXX
)

target_link_libraries(
  lattice_symmetries_c_sources
  PUBLIC 
    OpenMP::OpenMP_C
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
set(LATTICE_SYMMETRIES_PC_LIBDIR ${CMAKE_INSTALL_FULL_LIBDIR})
set(LATTICE_SYMMETRIES_PC_INCLUDEDIR ${CMAKE_INSTALL_FULL_INCLUDEDIR})
set(LATTICE_SYMMETRIES_PC_VERSION ${PROJECT_VERSION})
set(LATTICE_SYMMETRIES_PC_LIBS "-L\${libdir} -llattice_symmetries")
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
    "${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig"
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

if(${PROJECT_NAME}_ENABLE_PROFILING)
  add_subdirectory(profile)
endif()

#
# Benchmarks
#
# add_subdirectory(benchmark)
