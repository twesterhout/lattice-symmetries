cmake_minimum_required(VERSION 3.15)

project(${CMAKE_PROJECT_NAME}Profile LANGUAGES C CXX)

add_executable(${PROJECT_NAME}
    main.c
)

target_link_libraries(
  ${PROJECT_NAME} 
  PUBLIC
    lattice_symmetries
)
