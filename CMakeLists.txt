cmake_minimum_required(VERSION 3.16)
project(point_charge_symmetry)
include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://github.com/molpro/library-manager.git
        GIT_TAG 0.11.5
)
FetchContent_MakeAvailable(LibraryManager)

LibraryManager_Project()
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(dependencies)
add_subdirectory(src)

if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    enable_testing()
    add_subdirectory(test)
    add_subdirectory(doc)
endif ()