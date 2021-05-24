#include "Molecule.h"
#include <numeric>
#include <sstream>
#include <cmath>
#define N_PERIODIC_TABLE 105
const char* const periodic_table[N_PERIODIC_TABLE] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db"};

namespace molpro::point_charge_symmetry {
static std::string upper_string(std::string s) {
  transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return toupper(c); });
  return s;
}
Molecule::Molecule(const std::string& filename) {
  std::ifstream f(filename);
  int n;
  std::string title;
  std::getline(f, title);
  n = std::stoi(title);
  std::getline(f, m_title);
  std::string el;
  for (int i = 0; i < n; i++) {
    double x, y, z, q = 0;
    f >> el >> x >> y >> z;
    for (int e = 0; e < N_PERIODIC_TABLE; e++)
      if (upper_string(el) == upper_string(periodic_table[e]))
        q = e + 1;
    if (q == 0)
      q = std::stoi(el);
    m_atoms.emplace_back(Eigen::Vector3d{x, y, z}, q, el);
  }
}
Molecule::Molecule(const Eigen::MatrixXd& coordinates, const Eigen::VectorXd& charges) {
 for (int i=0; i<coordinates.cols(); i++) {
   m_atoms.emplace_back(coordinates.col(i), charges(i), std::to_string(std::lround(charges(i))));
 }
}

void Molecule::write(const std::string& filename, const std::string& title, const std::string& format) {
  std::ofstream f(filename);
//  Eigen::IOFormat fmt(10, 0);
  f.precision(10);
  f << std::fixed;
  f << m_atoms.size() << "\n" << (title.empty() ? m_title : title) << "\n";
  for (const auto& atom : m_atoms) {

//    f << atom.name << " " << atom.position.transpose().format(fmt) << std::endl;
  f << atom.name; for (int i=0; i<3; i++) f << " " << atom.position(i); f << std::endl;
  }
}

Eigen::Vector3d Molecule::centre_of_charge() const {
  Eigen::Vector3d result{0, 0, 0};
  for (const auto& atom : m_atoms)
    result += atom.position * atom.charge;
  result /= std::accumulate(m_atoms.begin(), m_atoms.end(), 0, [](auto init, auto atom) { return init + atom.charge; });
  return result;
}

Eigen::Matrix3d Molecule::inertia_tensor() const {
  Eigen::Matrix3d result = Eigen::Matrix3d::Zero();
  auto coc = centre_of_charge();
  for (const auto& atom : m_atoms) {
    result += atom.charge * ((atom.position - coc).dot(atom.position - coc) * Eigen::Matrix3d::Identity() -
                             (atom.position - coc) * (atom.position - coc).transpose());
    //    std::cout << "shifted atom position " << (atom.position - coc).transpose() << std::endl;
  }
  //  std::cout << "Inertia tensor\n" << result << std::endl;
  return result;
}

Eigen::Matrix3d Molecule::inertial_axes() const {
  Eigen::Matrix3d result;
  auto it = inertia_tensor();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(it);
  Eigen::Matrix3d ev = solver.eigenvectors().eval();
  for (int i = 0; i < 3; i++)
    if (ev(i, i) < 0)
      ev.col(i) = -ev.col(i);
  if (ev.determinant() < 0) {
    Eigen::Vector3d ev3 = ev.col(2).eval();
    ev.col(2) = ev.col(1);
    ev.col(1) = ev3;
  }
  //  std::cout << "Inertia tensor eigenvalues " << solver.eigenvalues().transpose() << std::endl;
  //  std::cout << "Inertia tensor eigenvectors\n" << ev << std::endl;
  return ev;
}

std::string Molecule::str() const {
  std::stringstream ss;
  ss << m_title << std::endl;
  for (const auto& atom : m_atoms)
    ss << atom.name << "(" << atom.charge << ") " << atom.position.transpose() << std::endl;
  ss << "Centre of charge: " << centre_of_charge().transpose() << std::endl;
  return ss.str();
}
} // namespace molpro::point_charge_symmetry
