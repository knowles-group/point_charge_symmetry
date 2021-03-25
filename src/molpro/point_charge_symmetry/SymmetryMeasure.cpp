#include "SymmetryMeasure.h"
#include <sstream>
namespace molpro::point_charge_symmetry {
Atom SymmetryMeasure::image(const Atom& source, const Operator& op) {
  return Atom(op(source.position), source.charge, source.name);
}
size_t SymmetryMeasure::image_neighbour(size_t atom_index, const Operator& op) {
  size_t result;
  double closest_distance = 1e50;
  //  std::cout << "Atom "<<m_molecule.m_atoms[atom_index].position <<std::endl;
  auto image = this->image(m_molecule.m_atoms[atom_index], op);
  //  std::cout << "Image "<<image.position <<std::endl;
  for (size_t i = 0; i < m_molecule.m_atoms.size(); i++) {
    if (m_molecule.m_atoms[i].charge == m_molecule.m_atoms[atom_index].charge) {
      double dist = (m_molecule.m_atoms[i].position - image.position).norm();
      if (dist < closest_distance) {
        closest_distance = dist;
        result = i;
      }
    }
  }
  return result;
}

double SymmetryMeasure::operator()(int operator_index) const {
  double result = 0;
  auto start = operator_index < 0 ? m_group.begin() : m_group.begin() + operator_index;
  auto end = operator_index < 0 ? m_group.end() : m_group.begin() + operator_index + 1;
  for (auto op = start; op < end; op++) {
//    std::cout << "Operator " << (*op)->name() << std::endl;
//    std::cout << "Frame parameters";
//    for (int i=0; i<6;i++)std::cout<<" "<<(*op)->coordinate_system().data()[i];
//    std::cout << std::endl;
    for (int a = 0; a < m_molecule.m_atoms.size(); a++) {
//            std::cout << "Atom " << a <<": "<< m_molecule.m_atoms[a].position.transpose() << std::endl;
//            std::cout << "Mapped Atom " << a <<": "<< (**op)(m_molecule.m_atoms[a].position).transpose() << std::endl;
      int ai = m_neighbours[op - m_group.begin()][a];
      //      std::cout << "Image " << ai << m_molecule.m_atoms[ai].position.transpose() << std::endl;
      auto dist = ((**op)(m_molecule.m_atoms[a].position) - m_molecule.m_atoms[ai].position).norm();
      //      std::cout << "Atom a dist=" << dist << std::endl;
      const auto zr = m_molecule.m_atoms[a].charge * dist;
      //      std::cout << "measure "<< 1 - std::exp(-zr) * (1 + zr + zr * zr / 3)<<std::endl;
      result += 1 - std::exp(-zr) * (1 + zr + zr * zr / 3);
    }
  }
//  std::cout << "result "<<result<<std::endl;
  return result;
}

CoordinateSystem::parameters_t SymmetryMeasure::coordinate_system_gradient(int operator_index) const {
  CoordinateSystem::parameters_t result{0, 0, 0, 0, 0, 0};
  const auto origin = m_group.coordinate_system().origin();
  const auto axes = m_group.coordinate_system().axes();
  const auto axes_gradient = m_group.coordinate_system().axes_gradient();
  auto start = operator_index < 0 ? m_group.begin() : m_group.begin() + operator_index;
  auto end = operator_index < 0 ? m_group.end() : m_group.begin() + operator_index + 1;
  for (auto op = start; op < end; op++) {
    //    std::cout << "Operator " << (*op)->name() << std::endl;
    for (int a = 0; a < m_molecule.m_atoms.size(); a++) {
      CoordinateSystem::parameters_t grad_d{0, 0, 0, 0, 0, 0};
//            std::cout << "Atom " << a << m_molecule.m_atoms[a].position.transpose() << std::endl;
//            std::cout << "Mapped Atom " << a << (**op)(m_molecule.m_atoms[a].position).transpose() << std::endl;
      int ai = m_neighbours[op - m_group.begin()][a];
//            std::cout << "Image " << ai << m_molecule.m_atoms[ai].position.transpose() << std::endl;
      auto d = ((**op)(m_molecule.m_atoms[a].position) - m_molecule.m_atoms[ai].position).eval();
//      std::cout << "d "<<d.transpose()<<std::endl;
      auto dist = d.norm();
//            std::cout << "Atom a dist=" << dist << std::endl;
      auto opgrad = (**op).operator_gradient(m_molecule.m_atoms[a].position,2,1e-4); // TODO analytic instead
//      std::cout << "d "<<d.transpose()<<std::endl;
        if (dist > 0)
        for (int i = 0; i < 6; i++) {
//          std::cout << "opgrad\n"<<opgrad[i]<<std::endl;
          for (int j = 0; j < 3; j++)
            grad_d[i] += opgrad[i][j] * d[j] / dist;
        }
      const auto zr = m_molecule.m_atoms[a].charge * dist;
      //      std::cout << "measure "<< 1 - std::exp(-zr) * (1 + zr + zr * zr / 3)<<std::endl;
      auto factor = std::exp(-zr) * m_molecule.m_atoms[a].charge * zr * (1 + zr) / 3;
      std::transform(grad_d.begin(), grad_d.end(), result.begin(), result.begin(),
                     [factor](const double& a, const double& b) { return b + a * factor; });
    }
  }
  return result;
}

std::string SymmetryMeasure::str() const {
  std::stringstream ss;
  ss << "SymmetryMeasure Group=" << m_group.name() << "\n";
  auto op = m_group.begin();
  for (const auto& nop : this->m_neighbours) {
    ss << "Image neighbours for symmetry operator " << (*op++)->name() << ":";
    for (const auto& image : nop)
      ss << " " << image;
    ss << std::endl;
  }
  return ss.str();
}

void SymmetryMeasure::optimise_frame() { std::cout << "optimise_frame" << std::endl; }

Group discover_group(const Molecule& molecule, double threshold) {
  Group result;
  return result;
}

} // namespace molpro::point_charge_symmetry
