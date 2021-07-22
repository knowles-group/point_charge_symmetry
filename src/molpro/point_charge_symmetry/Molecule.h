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
  /*!
   * @brief Add noise to the geometry, with each coordinate displaced randomly
   * @param amplitude Range of noise. Coordinate displacements are a random number in [-amplitude,amplitude)
   */
  void randomise(double amplitude);
};

inline std::ostream& operator<<(std::ostream& os, const Molecule& op) {
  os << op.str();
  return os;
}

double cartesian_distance(const Molecule& molecule1, const Molecule& molecule2);

/*!
 * @brief Compute a cooordinate-frame-independent measure of difference between two molecules.
 *
 * For each pair of atoms, we calculate for each molecule f0(rij) = 1-exp(-rij)(1+rij(1+rij/3))), which is rij*rij/6
 * for small rij, but tends asymptotically to 1 for large rij. This has the effect of avoiding overweighting of deviations in contributions between distant atoms.
 * The measure is then sum_{i>j} (f0(rij,A)-f0(rij,B))^2 / (n*(n-1)/2), which is bounded by 0,1
 * @param molecule
 * @param molecule0
 * @return
 */
std::pair<double, std::vector<double>> distance(const Molecule& molecule, const Molecule& molecule0);

} // namespace molpro::point_charge_symmetry

#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_MOLECULE_H_
