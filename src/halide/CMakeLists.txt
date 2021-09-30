# Generator
add_executable(generator generator.cpp)
target_link_libraries(generator PRIVATE Halide::Generator)

function(generate_kernels_for_arch)
  set(_options)
  set(_oneValueArgs
    FROM
    GENERATOR
    BASE_NAME
    BASE_ARCH
    OUTPUT_LIST
  )
  set(_multiValueArgs
    SUFFIXES
    FEATURES
    PARAMS
  )
  cmake_parse_arguments(generate_kernels_for_arch
      "${_options}" "${_oneValueArgs}" "${_multiValueArgs}" ${ARGN})

  # Base library
  set(_output_list)
  set(_base_library "${generate_kernels_for_arch_BASE_NAME}")
  add_halide_library(${_base_library}
      FROM ${generate_kernels_for_arch_FROM}
      GENERATOR ${generate_kernels_for_arch_GENERATOR}
      PARAMS ${generate_kernels_for_arch_PARAMS}
      TARGETS ${generate_kernels_for_arch_BASE_ARCH}
  )
  list(APPEND _output_list ${_base_library})

  set(_start 0)
  list(LENGTH generate_kernels_for_arch_SUFFIXES _stop)
  math(EXPR _stop "${_stop} - 1")
  if (${_stop} GREATER_EQUAL 0)
    foreach(_i RANGE ${_start} ${_stop})
      list(GET generate_kernels_for_arch_SUFFIXES ${_i} _suffix)
      list(GET generate_kernels_for_arch_FEATURES ${_i} _feature)
      set(_library "${generate_kernels_for_arch_BASE_NAME}${_suffix}")
      add_halide_library(${_library}
          FROM ${generate_kernels_for_arch_FROM}
          GENERATOR ${generate_kernels_for_arch_GENERATOR}
          PARAMS ${generate_kernels_for_arch_PARAMS}
          TARGETS ${generate_kernels_for_arch_BASE_ARCH}${_feature}
          USE_RUNTIME ${_base_library}
      )
      list(APPEND _output_list ${_library})
    endforeach()
  endif()
  set(${generate_kernels_for_arch_OUTPUT_LIST} ${_output_list} PARENT_SCOPE)
endfunction()


if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "(x86|AMD64)")
  set(base_arch x86-64-${os_for_Halide})
  set(arch_suffixes _sse41 _avx _avx2)
  set(arch_features -sse41 -avx-sse41 -avx-avx2-fma-sse41)
else()
  set(base_arch host)
  set(arch_suffixes)
  set(arch_features)
endif()

set(spin_inversion_values 0       1         -1)
set(kernel_names          general symmetric antisymmetric)
set(HalideKernels)

list(LENGTH spin_inversion_values _length)
math(EXPR _length "${_length} - 1")
foreach(_i RANGE ${_length})
  list(GET spin_inversion_values ${_i} _param)
  list(GET kernel_names ${_i} _name)

  generate_kernels_for_arch(
    FROM generator
    GENERATOR state_info_generator
    PARAMS spin_inversion=${_param}
    BASE_NAME ls_internal_state_info_${_name}_kernel_64
    BASE_ARCH ${base_arch}
    SUFFIXES ${arch_suffixes}
    FEATURES ${arch_features}
    OUTPUT_LIST _local_kernels
  )
  list(APPEND HalideKernels ${_local_kernels})

  generate_kernels_for_arch(
    FROM generator
    GENERATOR is_representative_generator
    PARAMS spin_inversion=${_param}
    BASE_NAME ls_internal_is_representative_${_name}_kernel_64
    BASE_ARCH ${base_arch}
    SUFFIXES ${arch_suffixes}
    FEATURES ${arch_features}
    OUTPUT_LIST _local_kernels
  )
  list(APPEND HalideKernels ${_local_kernels})
endforeach()
  
add_library(lattice_symmetries_halide_kernels STATIC kernels.hpp kernels.cpp)
target_link_libraries(lattice_symmetries_halide_kernels
  PRIVATE ${HalideKernels} Halide::Runtime
)
target_include_directories(lattice_symmetries_halide_kernels
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../src
)
target_include_directories(lattice_symmetries_halide_kernels
  SYSTEM PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party
)
target_compile_features(lattice_symmetries_halide_kernels PRIVATE cxx_std_17)
target_compile_definitions(lattice_symmetries_halide_kernels PRIVATE OUTCOME_DISABLE_EXECINFO)
set_project_warnings(lattice_symmetries_halide_kernels)
disable_rtti_and_exceptions(lattice_symmetries_halide_kernels)
set_property(TARGET lattice_symmetries_halide_kernels PROPERTY POSITION_INDEPENDENT_CODE ON)
