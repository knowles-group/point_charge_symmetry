include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://gitlab.com/molpro/LibraryManager.git
        GIT_TAG 0.7.2
)
FetchContent_MakeAvailable(LibraryManager)
LibraryManager_Project(point_charge_symmetry)
set(CMAKE_CXX_STANDARD 17)
LibraryManager_Add(point_charge_symmetry NAMESPACE molpro)

add_subdirectory(molpro)
LibraryManager_Install(point_charge_symmetry)

include(FetchContent)
FetchContent_Declare(tclap GIT_REPOSITORY https://gitlab.com/molpro/tclap.git)
FetchContent_GetProperties(tclap)
if (NOT tclap_POPULATED)
    FetchContent_Populate(tclap)
endif ()

add_executable(program program.cpp molpro/point_charge_symmetry/SymmetryMeasureC.cpp molpro/point_charge_symmetry/SymmetryMeasureC.h)
target_link_libraries(program PRIVATE molpro::point_charge_symmetry)
target_include_directories(program PRIVATE ${tclap_SOURCE_DIR}/include)
set_target_properties(program PROPERTIES OUTPUT_NAME point_charge_symmetry)
set_target_properties(program PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
set_target_properties(program PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
install(TARGETS program RUNTIME DESTINATION bin)