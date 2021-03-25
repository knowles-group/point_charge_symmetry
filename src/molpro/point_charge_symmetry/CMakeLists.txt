set_target_properties(point_charge_symmetry
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        )
LibraryManager_Append(point_charge_symmetry
        SOURCES Operator.cpp CoordinateSystem.cpp Molecule.cpp SymmetryMeasure.cpp
        PUBLIC_HEADER Operator.h CoordinateSystem.h Group.h Molecule.h SymmetryMeasure.h)
find_package(Eigen3 REQUIRED)
target_link_libraries(point_charge_symmetry PUBLIC Eigen3::Eigen)

include(FetchContent)
FetchContent_Declare(LinearAlgebra
        GIT_REPOSITORY https://gitlab.com/molpro/LinearAlgebra)
set(FORTRAN OFF)
set(MPI OFF)
set(LINEARALGEBRA_ARRAY_MPI3 OFF)
set(LINEARALGEBRA_ARRAY_GA OFF)
set(LINEARALGEBRA_ARRAY_HDF5 OFF)
FetchContent_MakeAvailable(LinearAlgebra)
target_link_libraries(point_charge_symmetry PUBLIC molpro::LinearAlgebra)