include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://gitlab+deploy-token-178345:RQZa2rmdwnV7Hdwg8XEq@gitlab.com/molpro/LibraryManager.git
        GIT_TAG master
)
FetchContent_MakeAvailable(LibraryManager)
LibraryManager_Project(point_charge_symmetry)
set(CMAKE_CXX_STANDARD 17)
LibraryManager_Add(point_charge_symmetry NAMESPACE molpro)

add_subdirectory(molpro)
LibraryManager_Install(point_charge_symmetry)