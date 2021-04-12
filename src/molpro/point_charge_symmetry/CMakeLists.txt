set_target_properties(point_charge_symmetry
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        )
LibraryManager_Append(point_charge_symmetry
        SOURCES Operator.cpp CoordinateSystem.cpp Molecule.cpp SymmetryMeasure.cpp SolverFactory-instantiation.cpp Group.cpp
        PUBLIC_HEADER Operator.h CoordinateSystem.h Group.h Molecule.h SymmetryMeasure.h)
find_package(Eigen3 REQUIRED)
target_link_libraries(point_charge_symmetry PUBLIC Eigen3::Eigen)

if (NOT TARGET molpro::LinearAlgebra)
    include(FetchContent)
    FetchContent_Declare(LinearAlgebra
            GIT_REPOSITORY https://gitlab.com/molpro/LinearAlgebra.git)
    set(FORTRAN OFF)
    set(MPI OFF)
    set(LINEARALGEBRA_ARRAY_MPI3 OFF)
    set(LINEARALGEBRA_ARRAY_GA OFF)
    set(LINEARALGEBRA_ARRAY_HDF5 OFF)
    FetchContent_MakeAvailable(LinearAlgebra)
endif ()
if (NOT TARGET molpro::Profiler)
    include(FetchContent)
    FetchContent_Declare(Profiler GIT_REPOSITORY https://gitlab.com/molpro/Profiler.git)
    FetchContent_MakeAvailable(Profiler)
endif ()
target_link_libraries(point_charge_symmetry PUBLIC molpro::LinearAlgebra molpro::Profiler)
