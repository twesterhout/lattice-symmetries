
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
    URL https://github.com/halide/Halide/releases/download/v12.0.1/Halide-12.0.1-${arch_for_Halide}-${bits_for_Halide}-${os_for_Halide}-5dabcaa9effca1067f907f6c8ea212f3d2b1d99a.${_archive_suffix}
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
