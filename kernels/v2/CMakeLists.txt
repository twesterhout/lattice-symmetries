cmake_minimum_required(VERSION 3.23)
project(lattice-symmetries-kernels VERSION 2.3 LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
  set(CMAKE_C_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES})
  set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

find_package(Halide REQUIRED)

# libffcall
find_library(LIBFFCALL_LIBRARY ffcall REQUIRED)
find_path(LIBFFCALL_INCLUDE_DIR ffcall-version.h REQUIRED)
add_library(libffcall INTERFACE)
target_link_libraries(libffcall INTERFACE ${LIBFFCALL_LIBRARY})
target_include_directories(libffcall INTERFACE ${LIBFFCALL_INCLUDE_DIR})

add_library(kernels_v2 STATIC kernels_v2.cpp reference.c indexing.c)
set_target_properties(kernels_v2 PROPERTIES PUBLIC_HEADER lattice_symmetries_types.h)
target_link_libraries(kernels_v2 PUBLIC Halide::Halide libffcall)
target_compile_features(kernels_v2 PRIVATE cxx_std_17)
target_compile_options(kernels_v2 PRIVATE
    # In case support for MSVC is ever added
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /WX>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic -Wno-cast-function-type>
)

add_executable(main main.cpp)
target_link_libraries(main PUBLIC kernels_v2)

include(GNUInstallDirs)
install(TARGETS kernels_v2 PUBLIC_HEADER)
