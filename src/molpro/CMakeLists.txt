include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://gitlab+deploy-token-178345:RQZa2rmdwnV7Hdwg8XEq@gitlab.com/molpro/LibraryManager.git
        GIT_TAG master
)
FetchContent_MakeAvailable(LibraryManager)
LibraryManager_Project(near_symmetry)
set(CMAKE_CXX_STANDARD 17)
LibraryManager_Add(near_symmetry NAMESPACE molpro)
LibraryManager_Append(near_symmetry SOURCES SymmetryOperator.cpp CoordinateSystem.cpp PUBLIC_HEADER SymmetryOperator.h CoordinateSystem.h)
find_package(Eigen3 REQUIRED)
target_link_libraries(near_symmetry PUBLIC Eigen3::Eigen)
add_library(molpro::near_symmetry ALIAS near_symmetry)