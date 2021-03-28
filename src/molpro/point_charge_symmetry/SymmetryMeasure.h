#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_SYMMETRYMEASURE_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_SYMMETRYMEASURE_H_
#include "Group.h"
#include "Molecule.h"

namespace molpro::point_charge_symmetry {

class SymmetryMeasure {
public:
  SymmetryMeasure(const Molecule& molecule, const Group& group) : m_molecule(molecule), m_group(group) {
    for (const auto& op : group) {
      m_neighbours.emplace_back();
      for (size_t i = 0; i < m_molecule.m_atoms.size(); i++)
        m_neighbours.back().push_back(image_neighbour(i, *op));
    }
  }

  double operator()(int operator_index = -1, int functional_form=0) const ;
  CoordinateSystem::parameters_t coordinate_system_gradient(int operator_index=-1, int functional_form=0) const;
  //  SymmetryMeasure(const Molecule& molecule, const Operator& op) : SymmetryMeasure(molecule, Group)
  std::string str() const;

  int optimise_frame();

protected:
  const Molecule& m_molecule;
  const Group& m_group;
  std::vector<std::vector<size_t>> m_neighbours;
  Atom image(const Atom& source, const Operator& op);
  size_t image_neighbour(size_t atom_index, const Operator& op);
};

inline std::ostream& operator<<(std::ostream& os, const SymmetryMeasure& sm) {
  os << sm.str();
  return os;
}

Group discover_group(const Molecule& molecule, double threshold=1e-10);

} // namespace molpro::point_charge_symmetry

#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_SYMMETRYMEASURE_H_
