cmake_minimum_required(VERSION 3.17)
project(point_charge_symmetry)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src)

include(FetchContent)
FetchContent_Declare(tclap GIT_REPOSITORY https://gitlab+deploy-token-103300:HKsBK8FZcZZshqAa_yyx@gitlab.com/molpro/tclap.git )
FetchContent_GetProperties(tclap)
if (NOT tclap_POPULATED)
    FetchContent_Populate(tclap)
endif ()
add_executable(symmetry_measure src/symmetry_measure.cpp)
target_link_libraries(symmetry_measure PRIVATE molpro::point_charge_symmetry)
target_include_directories(symmetry_measure PRIVATE ${tclap_SOURCE_DIR}/include)

if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    #    add_subdirectory(examples)
    enable_testing()
    add_subdirectory(test)
endif ()

# Documentation has to come last to make sure all dependencies have been added
#include(DependencyManagerDocs)

#DependencyManagerDocs_Add(${PROJECT_NAME}
#        FILES README.md
#        DOC_URL "https://molpro.gitlab.io/${PROJECT_NAME}/${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
#        )