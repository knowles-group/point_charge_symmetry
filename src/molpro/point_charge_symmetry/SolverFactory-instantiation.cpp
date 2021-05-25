#if defined(MOLPRO) && !defined(HAVE_MPI_H)
#define HAVE_MPI_H 1
#endif
#include "SymmetryMeasure.h"
#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
using Rvector = molpro::point_charge_symmetry::CoordinateSystem::parameters_t;
template class molpro::linalg::itsolv::SolverFactory<Rvector, Rvector>;
using Rvector2 = std::vector<double>;
template class molpro::linalg::itsolv::SolverFactory<Rvector2, Rvector2>;
