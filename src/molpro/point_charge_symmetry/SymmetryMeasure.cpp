#include "SymmetryMeasure.h"
#include <molpro/linalg/itsolv/OptimizeBFGS.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <sstream>
namespace molpro::point_charge_symmetry {
Atom SymmetryMeasure::image(const Atom& source, const Operator& op) {
  return Atom(op(source.position), source.charge, source.name);
}
size_t SymmetryMeasure::image_neighbour(size_t atom_index, const Operator& op) {
  size_t result;
  //  std::cout << "Operator "<<op<<std::endl;
  double closest_distance = 1e50;
  //    std::cout << "Atom "<<m_molecule.m_atoms[atom_index].position.transpose() <<std::endl;
  auto image = this->image(m_molecule.m_atoms[atom_index], op);
  //    std::cout << "Image "<<image.position.transpose() <<std::endl;
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

double SymmetryMeasure::operator()(int operator_index, int functional_form) const {
  double result = 0;
  auto start = operator_index < 0 ? m_group.begin() : m_group.begin() + operator_index;
  auto end = operator_index < 0 ? m_group.end() : m_group.begin() + operator_index + 1;
  for (auto op = start; op < end; op++) {
    //    std::cout << "Operator " << (*op)->name() << std::endl;
    //    std::cout << "Frame parameters";
    //    for (int i=0; i<6;i++)std::cout<<" "<<(*op)->coordinate_system().data()[i];
    //    std::cout << std::endl;
    for (int a = 0; a < m_molecule.m_atoms.size(); a++) {
      //                  std::cout << "Atom " << a <<": "<< m_molecule.m_atoms[a].position.transpose() << std::endl;
      //                  std::cout << "Mapped Atom " << a <<": "<< (**op)(m_molecule.m_atoms[a].position).transpose()
      //                  << std::endl;
      int ai = m_neighbours[op - m_group.begin()][a];
      //            std::cout << "Image " << ai << m_molecule.m_atoms[ai].position.transpose() << std::endl;
      auto dist = ((**op)(m_molecule.m_atoms[a].position) - m_molecule.m_atoms[ai].position).norm();
      //            std::cout << "Atom a dist=" << dist << std::endl;
      const auto zr = m_molecule.m_atoms[a].charge * dist;
      //            std::cout << "measure "<< 1 - std::exp(-zr) * (1 + zr + zr * zr / 3)<<std::endl;
      if (functional_form == 0)
        result += 1 - std::exp(-zr) * (1 + zr + zr * zr / 3);
      else if (functional_form == 1)
        result += zr * zr;
      else
        throw std::logic_error("invalid functional_form " + std::to_string(functional_form));
    }
  }
  //  std::cout << "result "<<result<<std::endl;
  return result;
}

CoordinateSystem::parameters_t SymmetryMeasure::coordinate_system_gradient(int operator_index,
                                                                           int functional_form) const {
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
      auto opgrad = (**op).operator_gradient(m_molecule.m_atoms[a].position, 2, 1e-4); // TODO analytic instead
      //      std::cout << "d "<<d.transpose()<<std::endl;
      if (dist > 0)
        for (int i = 0; i < 6; i++) {
          //          std::cout << "opgrad\n"<<opgrad[i]<<std::endl;
          for (int j = 0; j < 3; j++)
            grad_d[i] += opgrad[i][j] * d[j] / dist;
        }
      const auto zr = m_molecule.m_atoms[a].charge * dist;
      //      std::cout << "measure "<< 1 - std::exp(-zr) * (1 + zr + zr * zr / 3)<<std::endl;
      double factor;
      if (functional_form == 0)
        factor = std::exp(-zr) * m_molecule.m_atoms[a].charge * zr * (1 + zr) / 3;
      else if (functional_form == 1)
        factor = 2 * m_molecule.m_atoms[a].charge * zr;
      std::transform(grad_d.begin(), grad_d.end(), result.begin(), result.begin(),
                     [factor](const double& a, const double& b) { return b + a * factor; });
    }
  }
  constexpr bool numerical = false;
  if (numerical) {
    const double step = 1e-4;
    auto cs = m_group.coordinate_system();
    auto g = Group(cs, m_group);
    auto sm = SymmetryMeasure(m_molecule, g);
    CoordinateSystem::parameters_t resultn{0, 0, 0, 0, 0, 0};
    for (int i = 0; i < resultn.size(); i++) {
      cs.m_parameters[i] -= step;
      auto smm = sm(-1, functional_form);
      cs.m_parameters[i] += 2 * step;
      auto smp = sm(-1, functional_form);
      cs.m_parameters[i] -= step;
      resultn[i] = (smp - smm) / (2 * step);
    }
    std::cout << "*this " << *this << std::endl;
    std::cout << "sm " << sm << std::endl;
    std::cout << "resultn " << resultn << std::endl;
    std::cout << "result  " << result << std::endl;
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

int SymmetryMeasure::optimise_frame(CoordinateSystem& coordinate_system) {
  const double centre_of_charge_penalty = 0e0;
  const int centre_of_charge_penalty_power = 4;
  using Rvector = CoordinateSystem::parameters_t;
  //  std::cout << "optimise_frame" << std::endl;
  //  std::cout << "coordinate_system passed "&c
  //  coordinate_system.m_parameters={1.08791,0.778845,0.959933,0.823031,0.622547,0.611889};
  for (int c = 0; c < 0; c++) { // TODO remove testing only
    coordinate_system.m_parameters = {1, 1, 1, 1, 1, 1};
    double step = 1e-3;
    auto value0 = (*this)();
    auto grad0 = coordinate_system_gradient();
    coordinate_system.m_parameters[c] += step;
    auto valuep = (*this)();
    std::cout << "plus displacement parameters=";
    for (int i = 0; i < 6; i++)
      std::cout << " " << coordinate_system.m_parameters[i];
    std::cout << ", value=" << valuep << std::endl;
    coordinate_system.m_parameters[c] -= 2 * step;
    auto valuem = (*this)();
    std::cout << "minus displacement parameters=";
    for (int i = 0; i < 6; i++)
      std::cout << " " << coordinate_system.m_parameters[i];
    std::cout << ", value=" << valuem << std::endl;
    coordinate_system.m_parameters[c] += step;
    std::cout << "analytic=" << grad0[c] << ", numerical=" << (valuep - valuem) / (2 * step) << std::endl;
  }
  auto solver = molpro::linalg::itsolv::create_Optimize<Rvector, Rvector, Rvector>(
      "BFGS", "max_size_qspace=4,convergence_threshold=1e-8");
  int nwork = 1;
  for (int iter = 0; iter < 1000; iter++) {
    //    std::cout << coordinate_system << std::endl;
    auto value = (*this)(-1, 1);
    auto grad = coordinate_system_gradient(-1, 1);
    auto centre_of_charge_displacement = (coordinate_system.origin() - m_molecule.centre_of_charge()).eval();
    //    std::cout << grad[0]<<std::endl;
    //    std::cout << "coordinate_system.origin() " << coordinate_system.origin().transpose() << std::endl;
    //    std::cout << "m_molecule.centre_of_charge() " << m_molecule.centre_of_charge().transpose() << std::endl;
    //    std::cout << "centre_of_charge_displacement " << centre_of_charge_displacement.transpose() << std::endl;
    //        value=0; std::fill(grad.begin(),grad.end(),0);
    //    std::fill(grad.begin() + 3, grad.end(), 0);
    //    for (int i = 0; i < 3; i++)
    //      grad[i] = -grad[i];
    // std::fill(grad.begin(),grad.end()-3,0);
    //    std::cout << "gradient before penalty:";
    //    for (int i = 0; i < 6; i++)
    //      std::cout << " " << grad[i];
    //    std::cout << std::endl;

    //    value += centre_of_charge_penalty * std::pow(centre_of_charge_displacement.norm(),
    //    centre_of_charge_penalty_power); auto grad_penalty = centre_of_charge_penalty_power *
    //                        std::pow(centre_of_charge_displacement.norm(), centre_of_charge_penalty_power - 2) *
    //                        centre_of_charge_displacement;
    //    for (int i = 0; i < 3; i++)
    //      grad[i] += grad_penalty(i);

    //    std::cout << "gradient:";
    //    for (int i = 0; i < 6; i++)
    //      std::cout << " " << grad[i];
    //    std::cout << std::endl;
    //    std::cout << "before add_value "<<nwork;for (int i=0;i<6;i++)std::cout<<"
    //    "<<coordinate_system.m_parameters[i];std::cout<<std::endl; std::cout << "before add_value "<<nwork;for (int
    //    i=0;i<6;i++)std::cout<<" "<<grad[i];std::cout<<std::endl;
    std::cout << "value=" << value << std::endl;
    std::cout << "Current parameters " << coordinate_system.m_parameters << std::endl;
    auto current_parameters = coordinate_system.m_parameters;
    if (solver->add_value(coordinate_system.m_parameters, value, grad)) {
      std::cout << "after add_value " << nwork;
      for (int i = 0; i < 6; i++)
        std::cout << " " << grad[i];
      std::cout << std::endl;
      cout << "precondition" << std::endl;
      for (int i = 0; i < 6; i++)
        grad[i] /= 1;
    } else if (false){
      std::cout << "LINE SEARCH WAS SPECIFIED" << std::endl;
      auto new_parameters = coordinate_system.m_parameters;
      std::cout << "Original parameters " << current_parameters << std::endl;
      std::cout << "     New parameters " << new_parameters << std::endl;
      double dist = 0;
      for (int i = 0; i < 6; i++)
        dist += std::pow(current_parameters[i] - new_parameters[i], 2);
      dist = std::sqrt(dist);
      std::cout << "Distance " << dist << std::endl;
      std::cout << "SCAN" << std::endl;
      for (double x = -0.5; x < 1.1; x += .01) {
        for (int i = 0; i < 6; i++)
          coordinate_system.m_parameters[i] = (1 - x) * current_parameters[i] + x * new_parameters[i];
        int component = 5;
        auto fi = (*this)(-1, 1);
        coordinate_system.m_parameters[component] += 1e-2;
        auto fid = (*this)(-1, 1);
        coordinate_system.m_parameters[component] -= 1e-2;
        auto gradi = coordinate_system_gradient(-1, 1);
        std::cout << "component " << component << " numerical gradient=" << (fi - fid) / 1e-2
                  << ", analytical gradient=" << gradi[component] << std::endl;
        double grad_proj = 0;
        for (int i = 0; i < 6; i++)
          grad_proj += gradi[i] * (new_parameters[i] - current_parameters[i]);
        std::cout << "x=" << x << " value=" << (*this)(-1, 1) << ", projected gradient " << grad_proj << std::endl;
      }
      coordinate_system.m_parameters = new_parameters;
    }
    nwork = solver->end_iteration(coordinate_system.m_parameters, grad);
    std::cout << "after end_iteration " << nwork;
    for (int i = 0; i < 6; i++)
      std::cout << " " << coordinate_system.m_parameters[i];
    std::cout << std::endl;
    solver->report();
    if (nwork <= 0)
      return iter;
  }
  return -1;
}

Group discover_group(const Molecule& molecule, double threshold) {
  Group result;
  return result;
}

} // namespace molpro::point_charge_symmetry

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
using Rvector = molpro::point_charge_symmetry::CoordinateSystem::parameters_t;
template class molpro::linalg::itsolv::SolverFactory<Rvector, Rvector, Rvector>;
