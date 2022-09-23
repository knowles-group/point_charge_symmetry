set_target_properties(point_charge_symmetry
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        )
file(GENERATE OUTPUT molpro_config.h CONTENT "")
LibraryManager_Append(point_charge_symmetry
        SOURCES point_charge_symmetry_Operator.cpp CoordinateSystem.cpp Molecule.cpp SymmetryMeasure.cpp SolverFactory-instantiation.cpp Group.cpp Projector.cpp Projector.h SymmetryMeasureC.cpp
        PUBLIC_HEADER Operator.h CoordinateSystem.h Group.h Molecule.h SymmetryMeasure.h SymmetryMeasureC.h
        PRIVATE_HEADER ${CMAKE_CURRENT_BINARY_DIR}/molpro_config.h)
target_include_directories(point_charge_symmetry AFTER PUBLIC "${CMAKE_CURRENT_BINARY_DIR}")
LibraryManager_LAPACK(point_charge_symmetry PUBLIC Intel10_64lp_seq Intel10_64lp)
if (MKL_TYPE STREQUAL "lp64")
    message(STATUS "Eigen will use MKL")
    target_compile_definitions(point_charge_symmetry PUBLIC EIGEN_USE_MKL_ALL)
endif ()
find_package(Eigen3 REQUIRED)
target_link_libraries(point_charge_symmetry PUBLIC Eigen3::Eigen)

if (NOT TARGET molpro::iterative-solver)
    set(FORTRAN OFF)
    set(MPI OFF)
    set(LINEARALGEBRA_ARRAY_MPI3 OFF)
    set(LINEARALGEBRA_ARRAY_GA OFF)
    set(LINEARALGEBRA_ARRAY_HDF5 OFF)
    dependencymanager_Populate(iterative-solver)
endif ()
target_link_libraries(point_charge_symmetry PUBLIC molpro::iterative-solver)
