set_target_properties(point_charge_symmetry
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        )
LibraryManager_Append(point_charge_symmetry
        SOURCES SymmetryOperator.cpp CoordinateSystem.cpp
        PUBLIC_HEADER SymmetryOperator.h CoordinateSystem.h)
find_package(Eigen3 REQUIRED)
target_link_libraries(point_charge_symmetry PUBLIC Eigen3::Eigen)