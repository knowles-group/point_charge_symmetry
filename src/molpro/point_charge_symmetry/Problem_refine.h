#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_UTIL_PROBLEM_REFINE_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_UTIL_PROBLEM_REFINE_H_
#include "Projector.h"
#include "SymmetryMeasure.h"
#include <memory>
#include <molpro/linalg/itsolv/IterativeSolver.h>

namespace molpro::point_charge_symmetry {
class Problem_refine : public molpro::linalg::itsolv::Problem<std::vector<double>> {
  SymmetryMeasure& m_sm;
  Molecule& m_molecule;
  double m_distance_penalty = 0;
  std::unique_ptr<Projector> m_projector;
  const Molecule m_reference_molecule = Molecule{{}, {}};

public:
  Problem_refine(SymmetryMeasure& sm, Molecule& molecule, double distance_penalty = 0, bool project = false)
      : m_sm(sm), m_molecule(molecule), m_distance_penalty(std::move(distance_penalty)),
        m_projector(project ? std::make_unique<Projector>(m_sm.group(), m_molecule) : nullptr),
        m_reference_molecule(molecule) {
  }
  value_t residual(const container_t& parameters, container_t& residual) const override {
    size_t j = 0;
    for (auto& atom : m_molecule.m_atoms)
      for (int i = 0; i < 3; i++)
        atom.position(i) = parameters[j++];
    auto value = m_sm(-1, 1);
    residual = m_sm.atom_gradient(-1, 1);
    if (m_projector != nullptr) {
      std::cout << "parameters";
      for (const auto& e : parameters)
        std::cout << " " << e;
      std::cout << std::endl;
      std::cout << "residual before projection ";
      for (const auto& e : residual)
        std::cout << " " << e;
      std::cout << std::endl;
      m_projector->remove_symmetric(residual);
      std::cout << "residual after projection ";
      for (const auto& e : residual)
        std::cout << " " << e;
      std::cout << std::endl;
    }
    if (m_distance_penalty > 0) {
      auto dist = distance(m_molecule, m_reference_molecule);
      value += m_distance_penalty * dist.first;
      std::transform(residual.begin(), residual.end(), dist.second.begin(), residual.begin(),
                     [this](double a, double b) { return a + m_distance_penalty * b; });
    }
    return value;
  }
  bool diagonals(container_t& d) const override {
    std::fill(d.begin(), d.end(), 100);
    return true;
  }
};
} // namespace molpro::point_charge_symmetry

#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_UTIL_PROBLEM_REFINE_H_
