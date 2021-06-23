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
  Atom(Eigen::Vector3d r, double q, std::string name = "") : position(std::move(r)), charge(q), name(std::move(name)) {}
};

class Molecule {
public:
  std::vector<Atom> m_atoms;
  std::string m_title;
  explicit Molecule(const std::string& filename);
  explicit Molecule(const Eigen::MatrixXd& coordinates, const Eigen::VectorXd& charges);
  [[nodiscard]] std::string str() const;
  [[nodiscard]] Eigen::Vector3d centre_of_charge() const;
  [[nodiscard]] Eigen::Matrix3d inertia_tensor() const;
  [[nodiscard]] Eigen::Matrix3d inertial_axes() const;
  void write(const std::string& filename, const std::string& title = "", const std::string& format = "xyz");
  [[nodiscard]] size_t size() const { return m_atoms.size(); }
  [[nodiscard]] Eigen::Vector3d findaxis(int order) const;
};

inline std::ostream& operator<<(std::ostream& os, const Molecule& op) {
  os << op.str();
  return os;
}
} // namespace molpro::point_charge_symmetry

#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_MOLECULE_H_
