if(${PROJECT_NAME}_ENABLE_CLANG_TIDY)
    find_program(CLANGTIDY
      NAMES
        clang-tidy
        clang-tidy-11
        clang-tidy-10)
  if(CLANGTIDY)
    set(CMAKE_CXX_CLANG_TIDY ${CLANGTIDY} -extra-arg=-Wno-unknown-warning-option)
    message(STATUS "clang-tidy finished setting up.")
  else()
    message(SEND_ERROR "clang-cidy requested but executable not found.")
  endif()
endif()

if(${PROJECT_NAME}_ENABLE_CPPCHECK)
  find_program(CPPCHECK cppcheck)
  if(CPPCHECK)
    set(CMAKE_CXX_CPPCHECK ${CPPCHECK} --suppress=missingInclude --enable=all
                           --inline-suppr --inconclusive -i ${CMAKE_SOURCE_DIR}/imgui/lib)
    message(STATUS "Cppcheck finished setting up.")
  else()
    message(SEND_ERROR "Cppcheck requested but executable not found.")
  endif()
endif()
