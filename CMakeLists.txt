cmake_minimum_required(VERSION 3.17)
project(near_symmetry)

add_subdirectory(src)

add_executable(near_symmetry_main main.cpp)
target_link_libraries(near_symmetry_main PRIVATE molpro::near_symmetry)
