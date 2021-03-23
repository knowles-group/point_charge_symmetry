#include "Molecule.h"
#include <sstream>
#define N_PERIODIC_TABLE 105
const char *const periodic_table[N_PERIODIC_TABLE] = {
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
    "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
    "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db"
};

namespace molpro::point_charge_symmetry {
static std::string upper_string(std::string s) {
  transform(s.begin(), s.end(), s.begin(),
            [](unsigned char c) { return toupper(c); });
  return s;
}
Molecule::Molecule(const std::string &filename) {
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
      if (upper_string(el) == upper_string(periodic_table[e])) q = e + 1;
    m_atoms.emplace_back(Eigen::Vector3d{x, y, z}, q);
  }
}

std::string Molecule::str() const {
  std::stringstream ss;
  ss << m_title << std::endl;
  for (const auto &atom : m_atoms)
    ss << atom.q << " " << atom.r.transpose() << std::endl;
  return ss.str();
}
}
