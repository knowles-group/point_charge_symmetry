set_target_properties(point_charge_symmetry
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        )
LibraryManager_Append(point_charge_symmetry
        SOURCES point_charge_symmetry_Operator.cpp CoordinateSystem.cpp Molecule.cpp SymmetryMeasure.cpp SolverFactory-instantiation.cpp Group.cpp Projector.cpp Projector.h SymmetryMeasureC.cpp
        PUBLIC_HEADER Operator.h CoordinateSystem.h Group.h Molecule.h SymmetryMeasure.h SymmetryMeasureC.h)
LibraryManager_LAPACK(point_charge_symmetry PUBLIC Intel10_64lp_seq Intel10_64lp)
if (MKL_TYPE STREQUAL "lp64")
    message(STATUS "Eigen will use MKL")
    target_compile_definitions(${PROJECT_NAME} PUBLIC EIGEN_USE_MKL_ALL)
endif ()
find_package(Eigen3 REQUIRED)
target_link_libraries(point_charge_symmetry PUBLIC Eigen3::Eigen)

if (NOT TARGET molpro::LinearAlgebra)
#    include(FetchContent)
#    FetchContent_Declare(LinearAlgebra
#            GIT_REPOSITORY https://gitlab.com/molpro/LinearAlgebra.git GIT_TAG 34d83775a29ae864f3b9c64f15dfb4bbac9c8bf9)
    set(FORTRAN OFF)
    set(MPI OFF)
    set(LINEARALGEBRA_ARRAY_MPI3 OFF)
    set(LINEARALGEBRA_ARRAY_GA OFF)
    set(LINEARALGEBRA_ARRAY_HDF5 OFF)
#    FetchContent_MakeAvailable(LinearAlgebra)
    dependencymanager_Populate(LinearAlgebra)
endif ()
target_link_libraries(point_charge_symmetry PUBLIC molpro::LinearAlgebra)
#if (NOT TARGET molpro::Profiler)
#    include(FetchContent)
#    FetchContent_Declare(Profiler GIT_REPOSITORY https://gitlab.com/molpro/Profiler.git)
#    FetchContent_MakeAvailable(Profiler)
#endif ()
#target_link_libraries(point_charge_symmetry PUBLIC molpro::Profiler)
