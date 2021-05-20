#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_UTIL_PROBLEM_REFINE_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_UTIL_PROBLEM_REFINE_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/point_charge_symmetry/SymmetryMeasure.h>
namespace molpro::point_charge_symmetry {
class Problem_refine : public molpro::linalg::itsolv::Problem<std::vector<double>> {
  SymmetryMeasure& m_sm;
  Molecule& m_molecule;

public:
  Problem_refine(SymmetryMeasure& sm, Molecule& molecule) : m_sm(sm), m_molecule(molecule) {}
  value_t residual(const container_t& parameters, container_t& residual) const override {
    size_t j = 0;
    for (auto& atom : m_molecule.m_atoms)
      for (int i = 0; i < 3; i++)
        atom.position(i) = parameters[j++];
    residual = m_sm.atom_gradient(-1, 1);
    return m_sm(-1, 1);
  }
  bool diagonals(container_t& d) const override {
    std::fill(d.begin(), d.end(), 100);
    return true;
  }
};
}

#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_UTIL_PROBLEM_REFINE_H_
