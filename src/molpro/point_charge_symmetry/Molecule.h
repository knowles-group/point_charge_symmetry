#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_MOLECULE_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_MOLECULE_H_

#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
namespace molpro::point_charge_symmetry {
struct Atom {
  Eigen::Vector3d r;
  double q;
  Atom(Eigen::Vector3d r, double q) : r(r), q(q) {}
};
class Molecule {
 public:
  std::vector<Atom> m_atoms;
  std::string m_title;
  Molecule(const std::string &filename);
  std::string str() const;
};

inline std::ostream &operator<<(std::ostream &os, const Molecule &op) {
  os << op.str();
  return os;
}
}

#endif //POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_MOLECULE_H_
