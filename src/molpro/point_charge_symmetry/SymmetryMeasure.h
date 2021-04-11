#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_SYMMETRYMEASURE_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_SYMMETRYMEASURE_H_
#include "Group.h"
#include "Molecule.h"

namespace molpro::point_charge_symmetry {

class SymmetryMeasure {
public:
  SymmetryMeasure(const Molecule& molecule, const Group& group) : m_molecule(molecule), m_group(group) {
    reset_neighbours();
  }

  void reset_neighbours() {
    m_neighbours.clear();
    for (const auto& op : m_group) {
      m_neighbours.emplace_back();
      for (size_t i = 0; i < m_molecule.m_atoms.size(); i++)
        m_neighbours.back().push_back(image_neighbour(i, *op));
    }
  }

  double operator()(int operator_index = -1, int functional_form = 0, int verbosity = -1) const;
  [[nodiscard]] CoordinateSystem::parameters_t coordinate_system_gradient(int operator_index = -1, int functional_form = 0) const;
  [[nodiscard]] std::vector<double> atom_gradient(int operator_index = -1, int functional_form = 0) const;
  //  SymmetryMeasure(const Molecule& molecule, const Operator& op) : SymmetryMeasure(molecule, Group)
  [[nodiscard]] std::string str() const;

  void adopt_inertial_axes();
  int optimise_frame();

  bool spherical_top() {
    constexpr double tol = 1e-3;
    double large = 0, small = 1e20;
    for (int i = 0; i < 3; i++) {
      large = std::max(large, m_inertia_principal_values(i));
      small = std::min(small, m_inertia_principal_values(i));
    }
    if (large == 0)
      return true;
    return (1 - small / large) < tol;
  }

  bool symmetric_top() {
    constexpr double tol = 1e-3;
    return (std::abs(m_inertia_principal_values(0) - m_inertia_principal_values(1)) < tol or
            std::abs(m_inertia_principal_values(1) - m_inertia_principal_values(2)) < tol or
            std::abs(m_inertia_principal_values(0) - m_inertia_principal_values(2)) < tol) and
           not spherical_top();
  }

  [[nodiscard]] Molecule refine(int repeat=1) const;
  [[nodiscard]] CoordinateSystem::vec inertia_principal_values() const { return m_inertia_principal_values; }

protected:
  const Molecule& m_molecule;
  const Group& m_group;
  std::vector<std::vector<size_t>> m_neighbours;
  CoordinateSystem::vec m_inertia_principal_values = {1, 2, 3};
  Atom image(const Atom& source, const Operator& op);
  size_t image_neighbour(size_t atom_index, const Operator& op);
};

inline std::ostream& operator<<(std::ostream& os, const SymmetryMeasure& sm) {
  os << sm.str();
  return os;
}

Group discover_group(const Molecule& molecule, CoordinateSystem& coordinate_system, double threshold = 1e-10,
                     int verbosity = -1);
Group discover_group(const Molecule& molecule, double threshold = 1e-10, int verbosity = -1);

Molecule molecule_localised(const CoordinateSystem& coordinate_system, const Molecule& source);

} // namespace molpro::point_charge_symmetry

#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_SYMMETRYMEASURE_H_
