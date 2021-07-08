#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_SYMMETRYMEASUREC_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_SYMMETRYMEASUREC_H_

#include <cstddef>
/*!
 * @brief Evaluate the symmetry measure for given coordinates and group
 * @param atoms The number of point charges
 * @param groupname Schoenflies symbol for the point group
 * @param coordinates cartesian coordinates running fastest, ie Fortran dimension (3,atoms)
 * @param charges Charges of atoms, dimension (atoms)
 * @return The numerical value of the symmetry measure.
 */
extern "C" double SymmetryMeasureValue(const char* groupname, size_t atoms, const double* coordinates,
                                       const double* charges);

/*!
 * @brief Optimise the coordinate frame to minimise the symmetry measure for given coordinates and group
 * @param atoms The number of point charges
 * @param groupname Schoenflies symbol for the point group
 * @param coordinates cartesian coordinates running fastest, ie Fortran dimension (3,atoms). On exit they are
 * changed to values expressing position in the new coordinate frame, but relative distances and angles are
 * preserved.
 * @param charges Charges of atoms, dimension (atoms)
 * @return 0 if successful
 */
extern "C" int SymmetryMeasureOptimiseFrame(const char* groupname, size_t atoms, double* coordinates,
                                            const double* charges);

/*!
 * @brief Discover the highest order point group for a set of charges
 * @param Threshold for symmetry measure in order to accept membership of group
 * @param atoms The number of point charges
 * @param coordinates cartesian coordinates running fastest, ie Fortran dimension (3,atoms). On exit they are
 * changed to values expressing position in the new coordinate frame, but relative distances and angles are
 * preserved.
 * @param charges Charges of atoms, dimension (atoms)
 * @return Schoenflies symbol for the point group
 */
extern "C" char* SymmetryMeasureDiscoverGroup(double threshold, size_t atoms, double* coordinates,
                                              const double* charges);

/*!
 * @brief Optimise the coordinates to minimise the symmetry measure for given group
 * @param atoms The number of point charges
 * @param groupname Schoenflies symbol for the point group
 * @param coordinates cartesian coordinates running fastest, ie Fortran dimension (3,atoms). On exit they are
 * changed to values expressing position in the new coordinate frame, and minimising the symmetry measure; relative
 * distances and angles are not preserved.
 * @param charges Charges of atoms, dimension (atoms)
 * @return The value of the symmetry measure after refinement
 */
extern "C" void SymmetryMeasureRefine(const char* groupname, size_t atoms, double* coordinates, const double* charges);

#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_SYMMETRYMEASUREC_H_
