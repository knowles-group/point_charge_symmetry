cmake_minimum_required(VERSION 3.17)
project(point_charge_symmetry)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src)

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