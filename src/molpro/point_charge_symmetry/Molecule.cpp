#include "Molecule.h"
#include "Problem_refine.h"
#include "SymmetryMeasure.h"
#include <cmath>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <numeric>
#include <random>
#include <regex>
#include <sstream>
#include <unsupported/Eigen/MatrixFunctions>
#define N_PERIODIC_TABLE 106
const char* const periodic_table[N_PERIODIC_TABLE] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl",
    "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",
    "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db"};

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
    double x, y, z, q = -1;
    f >> el >> x >> y >> z;
    for (int e = 0; e < N_PERIODIC_TABLE; e++)
      if (upper_string(el) == upper_string(periodic_table[e]))
        q = e;
    if (q == -1)
      q = std::stoi(el);
    if (q == 0)
      q = 0.99; // charge of zero does not play well
    m_atoms.emplace_back(Eigen::Vector3d{x, y, z}, q, el);
  }
}
Molecule::Molecule(const Eigen::MatrixXd& coordinates, const Eigen::VectorXd& charges) {
  for (int i = 0; i < coordinates.cols(); i++) {
    m_atoms.emplace_back(coordinates.col(i), charges(i) == 0 ? 0.99 : charges(i), // charge of zero does not play well
                         std::to_string(std::lround(charges(i))));
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
    f << atom.name;
    for (int i = 0; i < 3; i++)
      f << " " << atom.position(i);
    f << std::endl;
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
static bool coplanar(const Eigen::Vector3d p1, const Eigen::Vector3d p2, const Eigen::Vector3d p3,
                     const Eigen::Vector3d p4, double threshold = 1e-6) {
  Eigen::Matrix3d discriminant;
  discriminant.col(0) = p2 - p1;
  discriminant.col(1) = p3 - p1;
  discriminant.col(2) = p4 - p1;
  //    std::cout << "planarity discriminant "<<discriminant.determinant()<<std::endl;
  return std::abs(discriminant.determinant()) < threshold;
}
Eigen::Vector3d Molecule::findaxis(int order) const {
  Eigen::Vector3d result;
  std::vector<bool> distant(this->size());
  std::vector<double> distance;
  auto centre = this->centre_of_charge();
  const double distance_threshold = 1e-3;
  for (const auto& atom : this->m_atoms) {
    const auto dist = (atom.position - centre).norm();
    distance.push_back(dist > distance_threshold ? dist : std::numeric_limits<double>::max());
  }
  for (int i = 0; i < 3; ++i) { // first of all try the raw coordinate axes
    CoordinateSystem cs;
    result = Eigen::Matrix3d::Identity().col(i);
    auto group = Group(cs);
    group.add(Rotation(group.coordinate_system(), result, order));
    SymmetryMeasure sm(*this, group);
    if (sm() < 1e-3)
      return result;
  }
  std::vector<size_t> atoms;
  size_t atom1 = std::min_element(distance.begin(), distance.end()) - distance.begin(); // the nearest distant atom
  { // next see if this atom could itself define the axis
    result = this->m_atoms[atom1].position - centre;
    CoordinateSystem cs;
    auto group = Group(cs);
    group.add(Rotation(group.coordinate_system(), result, order));
    SymmetryMeasure sm(*this, group);
    //    std::cout << "Trying nearest as axis definer " << result.transpose() << ", symmetry measure = " << sm()
    //              << std::endl;
    if (sm() < 1e-3)
      return result;
  }
  atoms.push_back(atom1);
  //    std::cout << "findaxis atom " << atom1 << " " << m_atoms[atom1].position.transpose() << std::endl;
  using neighbour_t = std::pair<size_t, double>;
  std::set<neighbour_t> near_neighbours;
  neighbour_t nearest_neighbour{0, std::numeric_limits<double>::max()};
  for (size_t candidate = 0; candidate < this->size(); candidate++)
    if (this->m_atoms[candidate].charge == this->m_atoms[atom1].charge and candidate != atom1) {
      const auto dist = (this->m_atoms[candidate].position - this->m_atoms[atom1].position).norm();
      near_neighbours.insert({candidate, dist});
      if (dist < nearest_neighbour.second)
        nearest_neighbour = {candidate, dist};
    }
  for (auto neighbour = near_neighbours.begin(); neighbour != near_neighbours.end();)
    if (neighbour->second > nearest_neighbour.second * 2.01)
      neighbour = near_neighbours.erase(neighbour);
    else
      neighbour++;

  //  std::cout << "nearest_neighbour " << nearest_neighbour.first << " : " << nearest_neighbour.second << std::endl;
  //  for (const auto& n : near_neighbours)
  //    std::cout << "near_neighbour " << n.first << " : " << n.second << std::endl;
  // find the expected regular polygon angle subtended at atom1 by a pair of neighbours
  size_t atom2 = atom1, atom3 = atom1;
  const auto pi = std::acos(double(-1));
  const auto expected_angle = pi - 2 * pi / order;
  const auto angle_tolerance = expected_angle / 100;
  for (const auto& n1 : near_neighbours)
    for (const auto& n2 : near_neighbours)
      if (n1 != n2) {
        const auto& n1xyz = this->m_atoms[n1.first].position;
        const auto& n2xyz = this->m_atoms[n2.first].position;
        const auto& n1vec = n1xyz - this->m_atoms[atom1].position;
        const auto& n2vec = n2xyz - this->m_atoms[atom1].position;
        const auto& angle = std::acos(n1vec.dot(n2vec) / (n1vec.norm() * n2vec.norm()));
        //        std::cout << n1.first << n2.first << " angle=" << angle << ", expected_angle=" << expected_angle <<
        //        std::endl;
        if (std::abs(angle - expected_angle) < angle_tolerance) {
          atom2 = n1.first;
          atom3 = n2.first;
        }
      }
  //    std::cout << "atom1, atom2, atom3 " << atom1 << " " << atom2 << " " << atom3 << std::endl;
  if (atom2 == atom1)
    return {0, 0, 0};
  atoms.push_back(atom2);
  atoms.push_back(atom3);

  for (int n = 3; n < order; n++) {
    //    std::cout << "Place polygon atom # "<<n<<std::endl;
    size_t atomn = 99; // the next atom closest to atom1 and of the same type, and within tolerance of the same distance
                       // from the origin
    //    std::cout << "n="<<n<<" atoms so far =";for (const auto& atom : atoms) std::cout << " "<<atom; std::cout <<
    //    std::endl ;
    for (size_t candidate = 0; candidate < this->size(); candidate++) {
      //      std::cout << "trying initial candidate "<<candidate<<std::endl;
      if (std::find(atoms.begin(), atoms.end(), candidate) == atoms.end() and
          this->m_atoms[candidate].charge == this->m_atoms[atom1].charge and
          coplanar(this->m_atoms[atom1].position, this->m_atoms[atoms[std::min(atoms.size() - 1, size_t(1))]].position,
                   this->m_atoms[atoms[std::min(atoms.size() - 1, size_t(2))]].position,
                   this->m_atoms[candidate].position, 1e-2)) {
        if (n == 1) {
          // for first neighbour, need to ensure there's a second neighbour at a similar distance
          // construct
        }
        atomn = candidate;
      }
    }
    //        std::cout << "setting initial atomn to "<<atomn<<std::endl;
    //        for (size_t candidate = 0; candidate < this->size(); candidate++)
    //          std::cout << "candidate "<<candidate<<" distance from 1 "<<(this->m_atoms[candidate].position -
    //          this->m_atoms[atom1].position).norm()<<" , distance from centre  "<<distance[candidate]<<", coplanar "<<
    //          (n<3?0:coplanar(this->m_atoms[atom1].position, this->m_atoms[atoms[std::min(atoms.size() - 1,
    //          size_t(1))]].position,
    //                                                                                                                                   this->m_atoms[atoms[std::min(atoms.size() - 1, size_t(2))]].position,
    //                                                                                                                                   this->m_atoms[candidate].position, 1e-2))<<std::endl;
    //    std::cout << "initial atomn="<<atomn<<std::endl;
    for (size_t candidate = 0; candidate < this->size(); candidate++)
      if (std::find(atoms.begin(), atoms.end(), candidate) == atoms.end() and
          distance[candidate] > distance[atom1] - 1e-2 and
          this->m_atoms[candidate].charge == this->m_atoms[atom1].charge and
          (this->m_atoms[candidate].position - this->m_atoms[atoms[n - 1]].position).norm() <
              (this->m_atoms[atomn].position - this->m_atoms[atoms[n - 1]].position).norm() + 1e-2 and
          coplanar(this->m_atoms[atom1].position, this->m_atoms[atoms[std::min(atoms.size() - 1, size_t(1))]].position,
                   this->m_atoms[atoms[std::min(atoms.size() - 1, size_t(2))]].position,
                   this->m_atoms[candidate].position, 1e-2) and
          (n < 2 or (this->m_atoms[candidate].position - this->m_atoms[atoms[n - 1]].position)
                            .dot(this->m_atoms[atoms[n - 1]].position - this->m_atoms[atoms[n - 2]].position) >
                        (this->m_atoms[atomn].position - this->m_atoms[atoms[n - 1]].position)
                            .dot(this->m_atoms[atoms[n - 1]].position - this->m_atoms[atoms[n - 2]].position)))
        atomn = candidate;
    atoms.push_back(atomn);
    //    std::cout << "findaxis atom " << atomn << " " << m_atoms[atomn].position.transpose() << std::endl;
  }
  result = (this->m_atoms[atoms[0]].position - this->m_atoms[atoms[1]].position)
               .cross(this->m_atoms[atoms[0]].position - this->m_atoms[atoms[2]].position);
  //    std::cout << "plane normal " << result.transpose() / result.norm() << std::endl;
  //    for (const auto& a : atoms)
  //      for (const auto& b : atoms)
  //    std::cout << "distance between atoms "<<a << ", "<<b<<" =
  //    "<<(m_atoms[a].position-m_atoms[b].position).norm()<<std::endl;
  return result / result.norm();
}
void Molecule::randomise(double amplitude) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(-amplitude, amplitude);
  for (auto& atom : m_atoms)
    for (int i = 0; i < 3; i++)
      atom.position(i) += dis(gen);
}

double cartesian_distance(const Molecule& molecule1, const Molecule& molecule2) {
  double dist = 0;
  for (size_t i = 0; i < molecule1.size(); ++i)
    dist += (molecule1.m_atoms[i].position - molecule2.m_atoms[i].position)
                .dot(molecule1.m_atoms[i].position - molecule2.m_atoms[i].position);
  dist = std::sqrt(dist);
  return dist;
}
std::pair<double, std::vector<double>> distance(const Molecule& molecule, const Molecule& molecule0) {
  std::vector<double> gradient(3 * molecule.size(), 0);
  double result = 0;
  auto distmod = [](auto& posi, auto& posj) {
    auto r = (posi - posj).norm();
    return (1 - std::exp(-r) * (1 + r * (1 + r / 3)));
  };
  auto distmodgrad = [](auto& posi, auto& posj) {
    auto r = (posi - posj).norm();
    return std::exp(-r) * r * (1 + r) / 3;
  };
  const auto scalefac = 1.0 / (molecule.size() * (molecule.size() - 1) / 2);
  for (size_t A = 0; A < molecule.size(); ++A) {
    for (size_t B = 0; B < A; ++B) {
      const auto posA = molecule.m_atoms[A].position;
      const auto posB = molecule.m_atoms[B].position;
      const auto pos0A = molecule0.m_atoms[A].position;
      const auto pos0B = molecule0.m_atoms[B].position;
      result += std::pow(distmod(posA, posB) - distmod(pos0A, pos0B), 2) * scalefac;
      //      std::cout << "interatomic distances "<<distmod(pos0A,pos0B)<<" -> "<<distmod(posA,posB)<<std::endl;
      for (size_t alpha = 0; alpha < 3; ++alpha)
        if ((posA - posB).norm() > 1e-10) {
          gradient[A * 3 + alpha] += 2 * scalefac * (distmod(posA, posB) - distmod(pos0A, pos0B)) *
                                     distmodgrad(posA, posB) * (posA(alpha) - posB(alpha)) / (posA - posB).norm();
          gradient[B * 3 + alpha] += 2 * scalefac * (distmod(posB, posA) - distmod(pos0B, pos0A)) *
                                     distmodgrad(posB, posA) * (posB(alpha) - posA(alpha)) / (posB - posA).norm();
        }
    }
  }
  //  std::cout << "result "<<result<<std::endl;
  return {result, gradient};
}

} // namespace molpro::point_charge_symmetry
