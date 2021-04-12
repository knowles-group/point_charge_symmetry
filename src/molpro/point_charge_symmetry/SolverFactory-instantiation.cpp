#include "SymmetryMeasure.h"
#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
using Rvector = molpro::point_charge_symmetry::CoordinateSystem::parameters_t;
template class molpro::linalg::itsolv::SolverFactory<Rvector, Rvector, Rvector>;
using Rvector2 = std::vector<double>;
template class molpro::linalg::itsolv::SolverFactory<Rvector2, Rvector2, Rvector2>;
