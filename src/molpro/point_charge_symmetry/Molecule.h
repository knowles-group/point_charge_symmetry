#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_MOLECULE_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_MOLECULE_H_

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <vector>
namespace molpro::point_charge_symmetry {
struct Atom {
  Eigen::Vector3d position;
  double charge;
  std::string name;
  Atom(Eigen::Vector3d r, double q, std::string name = "") : position(r), charge(q), name(std::move(name)) {}
};

class Molecule {
public:
  std::vector<Atom> m_atoms;
  std::string m_title;
  Molecule(const std::string& filename);
  std::string str() const;
};

inline std::ostream& operator<<(std::ostream& os, const Molecule& op) {
  os << op.str();
  return os;
}
} // namespace molpro::point_charge_symmetry

#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_MOLECULE_H_
