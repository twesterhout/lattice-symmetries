
# Detect operating system for Halide targets
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(os_for_Halide "linux")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(os_for_Halide "osx")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(os_for_Halide "windows")
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(bits_for_Halide 64)
else() 
    set(bits_for_Halide 32)
endif()

if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "(x86|AMD64)")
    set(arch_for_Halide x86)
else()
    set(arch_for_Halide arm)
endif()

# Find Halide
if(${PROJECT_NAME}_USE_SYSTEM_HALIDE)
    find_package(Halide REQUIRED COMPONENTS static)
else()
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
      set(_archive_suffix "zip")
  else()
      set(_archive_suffix "tar.gz")
  endif()
  FetchContent_Declare(
    HalideBinaryRelease
    URL https://github.com/halide/Halide/releases/download/v14.0.0/Halide-14.0.0-${arch_for_Halide}-${bits_for_Halide}-${os_for_Halide}-6b9ed2afd1d6d0badf04986602c943e287d44e46.${_archive_suffix}
    # URL https://github.com/halide/Halide/releases/download/v13.0.2/Halide-13.0.2-${arch_for_Halide}-${bits_for_Halide}-${os_for_Halide}-cf8c8f22eb507aedeba5a44f8ee20bc63757dc57.${_archive_suffix}
  )
  FetchContent_GetProperties(HalideBinaryRelease)
  if(NOT ${halidebinaryrelease_POPULATED})
    message(STATUS "[lattice-symmetries] Downloading binary release of Halide. This may take a while...")
    FetchContent_Populate(HalideBinaryRelease)
  endif()
  FetchContent_GetProperties(HalideBinaryRelease)
  find_package(ZLIB REQUIRED)
  find_package(HalideHelpers REQUIRED PATHS ${halidebinaryrelease_SOURCE_DIR}/lib/cmake)
  find_package(Halide REQUIRED COMPONENTS static PATHS ${halidebinaryrelease_SOURCE_DIR}/lib/cmake)
endif()
