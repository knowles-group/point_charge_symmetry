#if defined(MOLPRO) && !defined(HAVE_MPI_H)
#define HAVE_MPI_H 1
#endif
#include "SymmetryMeasure.h"
//#include <molpro/Profiler.h>
#include "Molecule.h"
#include "Problem_refine.h"
#include <cstdlib>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <regex>
#include <sstream>
#include <unsupported/Eigen/MatrixFunctions>
namespace molpro::point_charge_symmetry {
Atom SymmetryMeasure::image(const Atom& source, const Operator& op) const {
  return Atom(op(source.position), source.charge, source.name);
}
static inline size_t image_neighbour_inline(SymmetryMeasure& sm, size_t atom_index, const Operator& op) {
  //  auto prof = molpro::Profiler::single()->push("SymmetryMeasure::image_neighbour");
  const size_t no_result = 1000000;
  size_t result = no_result;
  //  std::cout << "Operator "<<op<<std::endl;
  double closest_distance = 1e50;
  //    std::cout << "Atom "<<m_molecule.m_atoms[atom_index].position.transpose() <<std::endl;
  auto image = sm.image(sm.m_molecule.m_atoms[atom_index], op);
  //    std::cout << "Image "<<image.position.transpose() <<std::endl;
  for (size_t i = 0; i < sm.m_molecule.m_atoms.size(); i++) {
    if (sm.m_molecule.m_atoms[i].charge == sm.m_molecule.m_atoms[atom_index].charge) {
      //      double dist = (sm.m_molecule.m_atoms[i].position - image.position).norm();
      //      double dist = (sm.m_molecule.m_atoms[i].position - image.position).dot((sm.m_molecule.m_atoms[i].position
      //      - image.position));
      double dist = std::pow(sm.m_molecule.m_atoms[i].position(0) - image.position(0), 2) +
                    std::pow(sm.m_molecule.m_atoms[i].position(1) - image.position(1), 2) +
                    std::pow(sm.m_molecule.m_atoms[i].position(2) - image.position(2), 2);
      if (dist < closest_distance) {
        closest_distance = dist;
        result = i;
      }
    }
  }
  if (result == no_result) {
    std::cout << "atom_index " << atom_index << std::endl;
    std::cout << "Atom " << sm.m_molecule.m_atoms[atom_index].position.transpose() << std::endl;
    std::cout << "op " << op << " " << op.str("Operator", true) << std::endl;
    std::cout << "Image " << image.position.transpose() << std::endl;
    for (size_t i = 0; i < sm.m_molecule.m_atoms.size(); i++) {
      std::cout << "atom " << sm.m_molecule.m_atoms[i].position.transpose() << std::endl;
    }
    throw std::logic_error("cannot get neighbour");
  }
  return result;
}
size_t SymmetryMeasure::image_neighbour(size_t atom_index, const Operator& op) {
  return image_neighbour_inline(*this, atom_index, op);
}

void SymmetryMeasure::reset_neighbours() {
  auto prof = molpro::Profiler::single()->push("SymmetryMeasure::reset_neighbours");
  m_neighbours.clear();
  size_t n = m_molecule.m_atoms.size();
  for (const auto& op : m_group) {
    m_neighbours.emplace_back();
    m_neighbours.back().reserve(n);
    for (size_t i = 0; i < n; i++)
      m_neighbours.back().emplace_back(image_neighbour_inline(*this, i, *op));
  }
}

double SymmetryMeasure::operator()(int operator_index, int functional_form, int verbosity) const {
  auto p = molpro::Profiler::single()->push("SymmetryMeasure()");
  double result = 0;
  auto start = operator_index < 0 ? m_group.begin() : m_group.find(m_group[operator_index]);
  auto end = operator_index < 0 ? m_group.end() : m_group.find(m_group[operator_index]);
  operator_index = std::max(operator_index, 0);
  if (verbosity > -1)
    std::cout << "!! SymmetryMeasure() " << std::endl;
  for (auto op = start; op != end; op++, ++operator_index) {
    if (verbosity > 0) {
      std::cout << "Operator " << (*op)->name() << std::endl;
      std::cout << "Frame parameters";
      for (int i = 0; i < 6; i++)
        std::cout << " " << (*op)->coordinate_system().data()[i];
      std::cout << std::endl;
    }
    for (size_t a = 0; a < m_molecule.m_atoms.size(); a++) {
      size_t ai = m_neighbours[operator_index][a];
      auto dist = ((**op)(m_molecule.m_atoms[a].position) - m_molecule.m_atoms[ai].position).norm();
      const auto zr = m_molecule.m_atoms[a].charge * dist;
      if (verbosity > 1) {

        std::cout << "Atom " << a << ": " << m_molecule.m_atoms[a].position.transpose() << std::endl;
        std::cout << "Mapped Atom " << a << ": " << (**op)(m_molecule.m_atoms[a].position).transpose() << std::endl;
        std::cout << "Image " << ai << m_molecule.m_atoms[ai].position.transpose() << std::endl;
        std::cout << "Atom a dist=" << dist << std::endl;
        std::cout << "measure " << 1 - std::exp(-zr) * (1 + zr + zr * zr / 3) << std::endl;
      }
      if (functional_form == 0)
        if (zr < 1e-8)
          result += zr * zr * (double(1) / 6 - zr * zr * (double(1) / 24 - zr / 45));
        else
          result += 1 - std::exp(-zr) * (1 + zr + zr * zr / 3);
      else if (functional_form == 1)
        result += zr * zr;
      else
        throw std::logic_error("invalid functional_form " + std::to_string(functional_form));
      if (verbosity > 0)
        std::cout << "new result " << result << std::endl;
    }
  }
  if (verbosity > -1)
    std::cout << "result " << result << std::endl;
  return result;
}

CoordinateSystem::parameters_t SymmetryMeasure::coordinate_system_gradient(int operator_index,
                                                                           int functional_form) const {
  auto p = molpro::Profiler::single()->push("SymmetryMeasure::coordinate_system_gradient()");
  CoordinateSystem::parameters_t result{0, 0, 0, 0, 0, 0};
  auto start = operator_index < 0 ? m_group.begin() : m_group.find(m_group[operator_index]);
  auto end = operator_index < 0 ? m_group.end() : m_group.find(m_group[operator_index]);
  operator_index = std::max(operator_index, 0);
  for (auto op = start; op != end; op++, ++operator_index) {
    //    std::cout << "Operator " << (*op)->name() << std::endl;
    for (size_t a = 0; a < m_molecule.m_atoms.size(); a++) {
      CoordinateSystem::parameters_t grad_d{0, 0, 0, 0, 0, 0};
      //            std::cout << "Atom " << a << m_molecule.m_atoms[a].position.transpose() << std::endl;
      //            std::cout << "Mapped Atom " << a << (**op)(m_molecule.m_atoms[a].position).transpose() << std::endl;
      size_t ai = m_neighbours[operator_index][a];
      //            std::cout << "Image " << ai << m_molecule.m_atoms[ai].position.transpose() << std::endl;
      auto d = ((**op)(m_molecule.m_atoms[a].position) - m_molecule.m_atoms[ai].position).eval();
      //      std::cout << "d "<<d.transpose()<<std::endl;
      auto dist = d.norm();
      //            std::cout << "Atom a dist=" << dist << std::endl;
      auto opgrad = (**op).operator_gradient(m_molecule.m_atoms[a].position, 0, 2e-3); // TODO analytic instead
      //      std::cout << "d "<<d.transpose()<<std::endl;
      if (dist > 0)
        for (int i = 0; i < 6; i++) {
          //                    std::cout << "opgrad\n"<<opgrad[i].transpose()<<std::endl;
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
      else
        throw std::logic_error("impossible functional form");
      std::transform(grad_d.begin(), grad_d.end(), result.begin(), result.begin(),
                     [factor](const double& a, const double& b) { return b + a * factor; });
    }
  }
  return result;
}

std::vector<double> SymmetryMeasure::atom_gradient(int operator_index, int functional_form) const {
  auto p = molpro::Profiler::single()->push("SymmetryMeasure::atom_gradient()");
  std::vector<double> result(m_molecule.m_atoms.size() * 3, 0);
  const auto axes = m_group.coordinate_system().axes();
  auto start = operator_index < 0 ? m_group.begin() : m_group.find(m_group[operator_index]);
  auto end = operator_index < 0 ? m_group.end() : m_group.find(m_group[operator_index]);
  operator_index = std::max(operator_index, 0);
  for (auto op = start; op != end; op++, ++operator_index) {
    //    std::cout << "Operator " << (*op)->name() << std::endl;
    CoordinateSystem::mat w; // T_t u^\dagger
    for (int j = 0; j < 3; j++)
      w.col(j) = (*op)->operator_local(axes.row(j));
    for (size_t a = 0; a < m_molecule.m_atoms.size(); a++) {
      //            std::cout << "Atom " << a << m_molecule.m_atoms[a].position.transpose() << std::endl;
      //            std::cout << "Mapped Atom " << a << (**op)(m_molecule.m_atoms[a].position).transpose() << std::endl;
      size_t ai = m_neighbours[operator_index][a];
      //            std::cout << "Image " << ai << m_molecule.m_atoms[ai].position.transpose() << std::endl;
      auto d = ((**op)(m_molecule.m_atoms[a].position) - m_molecule.m_atoms[ai].position).eval();
      //      std::cout << "d "<<d.transpose()<<std::endl;
      auto dist = d.norm();
      //            std::cout << "Atom a dist=" << dist << std::endl;
      if (dist > 0) {
        const auto zr = m_molecule.m_atoms[a].charge * dist;
        double factor;
        if (functional_form == 0)
          factor = std::exp(-zr) * m_molecule.m_atoms[a].charge * zr * (1 + zr) / 3;
        else if (functional_form == 1)
          factor = 2 * m_molecule.m_atoms[a].charge * zr;
        else
          throw std::logic_error("impossible functional form");
        factor /= dist;
        Eigen::Map<Eigen::Vector3d>(&result[3 * a], 3) += factor * w.transpose() * axes.transpose() * d;
        Eigen::Map<Eigen::Vector3d>(&result[3 * ai], 3) -= factor * d;
      }
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

void SymmetryMeasure::adopt_inertial_axes() {
  auto prof = molpro::Profiler::single()->push("SymmetryMeasure::adopt_inertial_axes");
  auto& coordinate_system = m_group.coordinate_system();
  auto& parameters = m_group.coordinate_system_parameters();
  //    std::cout << "adopt_inertial_axes initial, group="<<m_group.name()<<"\n" << coordinate_system.axes() <<
  //    std::endl;
  { // initialise to inertial axes in some orientation or other
    CoordinateSystem temporary_coordinate_system(coordinate_system.m_rotation_parameter_type,
                                                 m_molecule.centre_of_charge(), m_molecule.inertial_axes());
    std::copy(temporary_coordinate_system.m_parameters.begin(), temporary_coordinate_system.m_parameters.end(),
              parameters.begin());
  }
  //    std::cout << "adopt_inertial_axes first inertial\n" << coordinate_system.axes() << std::endl;
  //    std::cout << "centre of charge " << m_molecule.centre_of_charge() << std::endl;
  //    std::cout <<"Group: "<<m_group<<std::endl;
  int best_axis = 0;
  double best_axis_sm = 1e50;
  //  bool symmetric_top = false;
  for (int principal_axis = 0; principal_axis < 6; principal_axis++) {
    reset_neighbours();
    //            std::cout << "try axes, principal_axis="<<principal_axis <<"\n" << coordinate_system.axes() <<
    //            std::endl;
    //    auto sm = SymmetryMeasure(m_molecule, m_group);
    auto local_inertia_tensor =
        (coordinate_system.axes().transpose() * m_molecule.inertia_tensor() * coordinate_system.axes()).eval();
    //    std::cout << "local_inertia_tensor\n" << local_inertia_tensor << std::endl;
    double measure;
    m_inertia_principal_values = local_inertia_tensor.diagonal().eval();
    //    std::cout << "local_inertia_tensor "<<local_inertia_tensor<<std::endl;
    //    std::cout << "m_inertia_principal_values "<<m_inertia_principal_values<<std::endl;
    //    std::cout << m_inertia_principal_values(0)-m_inertia_principal_values(1)<<std::endl;
    //    std::cout << m_inertia_principal_values(0)-m_inertia_principal_values(2)<<std::endl;
    //    std::cout << m_inertia_principal_values(1)-m_inertia_principal_values(2)<<std::endl;
    if (
        //         std::abs(local_inertia_tensor(0, 0) - local_inertia_tensor(1, 1)) < 1e-5 *
        //         std::abs(local_inertia_tensor(0, 0) - local_inertia_tensor(2, 2))
        std::abs(local_inertia_tensor(0, 0) - local_inertia_tensor(1, 1)) <
        1e-2 * std::min(local_inertia_tensor(2, 2), double(1))
        //         std::abs(vals[0]-vals[1]) < 1e-2 or std::abs(vals[2]-vals[1]) < 1e-2
    ) {
      //      std::cout << local_inertia_tensor(0, 0) << " " << local_inertia_tensor(1, 1) << std::endl;
      //            std::cout << "Symmetric top with axes\n" << coordinate_system.axes() << std::endl;
      measure = 0;
      //      symmetric_top = true;
    } else
      measure = (*this)();
    //            std::cout << "Atomic coordinates in local frame\n" << std::endl;
    //        for (const auto &atom : m_molecule.m_atoms)
    //          std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
    CoordinateSystem::vec out_of_plane{0, 0, 0};
    for (const auto& atom : m_molecule.m_atoms) {
      //      auto coords= coordinate_system.to_local(atom.position);
      out_of_plane += coordinate_system.to_local(atom.position).cwiseAbs();
      //      std::cout << coords.transpose()<<std::endl;
      //      std::cout << coordinate_system.to_local(atom.position).cwiseAbs().transpose()<<std::endl;
      //            std::cout << out_of_plane.transpose()<<std::endl;
      //      std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
    }
    // TODO implement all recommended tie breakers, such as planar C2v molecules in yz plane
    auto tie_breaker = out_of_plane.dot(CoordinateSystem::vec{3, 2, 1}) * 1e-6 / out_of_plane.norm();
    //    std::cout << "tie_breaker="<<tie_breaker<<std::endl;
    measure += tie_breaker;
    //        std::cout << "principal_axis=" << principal_axis << " sm=" << measure << std::endl;
    //    std::cout << "out_of_plane "<<out_of_plane.transpose()<<std::endl;
    if (measure < best_axis_sm) {
      best_axis_sm = measure;
      best_axis = principal_axis;
    }
    //    std::cout << "axes before cycle\n" << coordinate_system.axes() << std::endl;
    coordinate_system.cycle_axes();
  }
  for (int principal_axis = 0; principal_axis < best_axis; principal_axis++)
    coordinate_system.cycle_axes();
  reset_neighbours();
  //    std::cout << "chosen initial coordinate system: " << best_axis << "\n" << coordinate_system << std::endl;
  //    std::cout << "Atomic coordinates in local frame\n" << std::endl;
  //    for (const auto& atom : m_molecule.m_atoms)
  //      std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
}

class Problem_optimise_frame : public molpro::linalg::itsolv::Problem<CoordinateSystem::parameters_t> {
  SymmetryMeasure& m_sm;

public:
  Problem_optimise_frame(SymmetryMeasure& sm) : m_sm(sm) {}
  value_t residual(const container_t& parameters, container_t& residual) const override {
    constexpr bool optimize_origin = true;
    //    std::cout << "parameters: "<<parameters<<std::endl;
    //    residual = m_sm.coordinate_system_gradient(-1, 1);
    //    if (std::inner_product(residual.begin(), residual.end(), residual.begin(), double(0)) > 1e-6) {
    m_sm.reset_neighbours();
    residual = m_sm.coordinate_system_gradient(-1, 1);
    //    }
    if (not optimize_origin)
      std::fill(residual.begin(), residual.begin() + 3, 0);
    //    std::cout << "residual: "<<residual<<std::endl;
    //    std::cout << "value: "<<m_sm(-1,1)<<std::endl;
    return m_sm(-1, 1);
  }
  bool diagonals(container_t& d) const override {
    std::fill(d.begin(), d.end(), 100);
    return true;
  }
};

int SymmetryMeasure::refine_frame(int verbosity) {
  auto prof = molpro::Profiler::single()->push("SymmetryMeasure::refine_frame");
  auto& parameters = m_group.coordinate_system_parameters();
  using Rvector = CoordinateSystem::parameters_t;
  auto solver =
      molpro::linalg::itsolv::create_Optimize<Rvector, Rvector>("BFGS", "max_size_qspace=9,convergence_threshold=1e-4");
  if (verbosity > 0) {
    std::cout << "initial";
    for (int i = 0; i < 6; i++)
      std::cout << " " << parameters[i];
    std::cout << std::endl;
    std::cout << m_group.coordinate_system().axes() << std::endl;
    std::cout << "Atomic coordinates in current local frame\n" << std::endl;
    for (const auto& atom : m_molecule.m_atoms)
      std::cout << atom.name << ": " << m_group.coordinate_system().to_local(atom.position).transpose() << std::endl;
  }
  auto problem = Problem_optimise_frame(*this);
  CoordinateSystem::parameters_t grad;
  solver->set_verbosity(linalg::itsolv::Verbosity::None);
  if (verbosity > 1)
    solver->set_verbosity(linalg::itsolv::Verbosity::Iteration);
  if (verbosity > 2)
    solver->set_verbosity(linalg::itsolv::Verbosity::Detailed);
  solver->set_max_iter(100);
  auto result = solver->solve(parameters, grad, problem, false);
  //    std::cout << "solve finds";
  //    std::cout << parameters;
  //    std::cout << std::endl;
  return result ? solver->statistics().iterations : -1;
}

template <typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) {
  for (const auto& e : v)
    s << " " << e;
  return s;
}

Molecule SymmetryMeasure::refine(double distance_penalty, bool project, int repeat) const {
  auto prof = molpro::Profiler::single()->push("SymmetryMeasure::refine");
  auto molecule = molecule_localised(m_group.coordinate_system(), this->m_molecule);
  //  std::cout << "refine initial molecule\n"<<molecule<<std::endl;
  auto group = Group(m_group.name());
  //  std::cout << "refine initial group\n"<<group<<std::endl;
  for (int attempt = 0; attempt < repeat; attempt++) {
    auto sm = SymmetryMeasure(molecule, group);
    using Rvector = std::vector<double>;
    Rvector parameters;
    for (const auto& atom : molecule.m_atoms)
      for (int i = 0; i < 3; i++)
        parameters.push_back(atom.position(i));
    const int verbosity = -1;
    auto solver = molpro::linalg::itsolv::create_Optimize<Rvector, Rvector>(
        "BFGS", "max_size_qspace=12,convergence_threshold=1e-7");
    solver->set_verbosity(linalg::itsolv::Verbosity::None);
    if (verbosity >= 1)
      solver->set_verbosity(linalg::itsolv::Verbosity::Iteration);
    auto problem = Problem_refine(sm, molecule, distance_penalty, project);
    auto grad = parameters;
    solver->solve(parameters, grad, problem);
    size_t j = 0;
    for (auto& atom : molecule.m_atoms)
      for (int i = 0; i < 3; i++)
        atom.position(i) = parameters[j++];
    //    std::cout << "in refine, from problem object final symmetry measure " << problem.residual(parameters, grad)
    //              << std::endl;
    //    std::cout << "in refine, final symmetry measure " << sm(-1) << std::endl;
    //    std::cout << "in refine, final symmetry measure " << sm(-1, 0) << std::endl;
    //    std::cout << "in refine, final symmetry measure " << sm(-1, 1) << std::endl;
  }
  return molecule;
}

bool test_group(const Molecule& molecule, const Group& group, double threshold, int verbosity) {
  auto prof = molpro::Profiler::single()->push("SymmetryMeasure::test_group");
  if (verbosity >= 0)
    std::cout << "test_group " << group.name() << std::endl;
  if (verbosity > 0)
    std::cout << "test_group " << group << std::endl;
  //  auto p = molpro::Profiler::single()->push("SymmetryMeasure::test_group(" + group.name() + ")");
  SymmetryMeasure sm0(molecule, group);
  sm0.adopt_inertial_axes();
  if (group.name() == "R3")
    return molecule.size() == 1;
  else if (molecule.size() == 1)
    return false;
  if ((group.name() == "Dinfh" or group.name() == "Cinfv") and
      std::min({sm0.inertia_principal_values()(0), sm0.inertia_principal_values()(1),
                sm0.inertia_principal_values()(2)}) > 1e-3) {
    //    std::cout << "abandoning testing linear, inertial principal values " <<
    //    sm.inertia_principal_values().transpose() << std::endl;
    return false;
  }
  if (verbosity > 0) {
    std::cout << "initial measure " << sm0() << std::endl;
    std::cout << group.coordinate_system().axes() << std::endl;
  }
  if (group.name().substr(0, 1) == "I") {
  }
  // scan
  if (sm0.spherical_top()) {
    //    std::cout << "looking at spherical top, group=" << group.name() << std::endl;
    //    std::cout << molecule << std::endl;
    group.coordinate_system().from_axes(find_axis_frame(molecule, group));
    //    std::cout << "find_axis_frame gives " << group.coordinate_system() << std::endl;
    //    if (true) {
    //      auto cs = group.coordinate_system();
    //      Group grot(cs);
    //      grot.add(group.highest_rotation());
    //      SymmetryMeasure smrot(molecule, grot);
    //      std::cout << "symmetry measure from rotation only " << smrot() << std::endl;
    //      std::cout << "grot.coordinate_system "<<grot.coordinate_system()<<std::endl;
    //      std::cout << "group.coordinate_system "<<group.coordinate_system()<<std::endl;
    //    }
    //    for (size_t i; i < group.size(); i++)
    //      std::cout << "Symmetry measure " << group[i].name() << " : " << sm(i) << std::endl;
    if (true) {
      CoordinateSystem cs = group.coordinate_system();
      Group grouprot(cs);
      for (int index = 0; group.highest_rotation(true, index).order() > 1; index++) {
        //        std::cout << "consider\n"<<group.highest_rotation(true,index)<<std::endl;
        if (index == 0 or
            std::abs(group.highest_rotation().axis().dot(group.highest_rotation(true, index).axis())) < 0.99) {
          grouprot.add(group.highest_rotation(true, index));
          if (index > 0)
            break;
        }
      }
      //      std::cout << "grouprot " << group.name() << "\n" << grouprot << std::endl;
      SymmetryMeasure sm(molecule, grouprot);
      const int nscan = 360;
      //      std::cout << "nscan=" << nscan << std::endl;
      double best_measure = 1e50;
      CoordinateSystem::parameters_t best_parameters;
      //      auto parameter_ranges = grouprot.coordinate_system().rotation_generator_ranges();
      auto axis = group.highest_rotation().axis();
      //      std::cout << "spin axis = " << axis.transpose() << std::endl;
      double angle = std::acos(double(-1)) * 2 / (nscan);
      Eigen::Matrix3d genx, geny, genz;
      genx << 0, 0, 0, 0, 0, 1, 0, -1, 0;
      geny << 0, 0, -1, 0, 0, 0, 1, 0, 0;
      genz << 0, 1, 0, -1, 0, 0, 0, 0, 0;
      Eigen::Matrix3d rot = (angle * (axis[0] * genx + axis[1] * geny + axis[2] * genz)).exp();
      //      std::cout << "rotation matrix\n" << rot << std::endl;
      for (int xscan = 0; xscan < nscan; xscan++) {
        grouprot.coordinate_system().from_axes(grouprot.coordinate_system().axes() * rot);
        sm.reset_neighbours();
        auto measure = sm();
        if (verbosity > 1)
          std::cout << "try " << measure << " ( current best=" << best_measure << ") "
                    << grouprot.coordinate_system_parameters() << std::endl;
        if (measure < best_measure) {
          best_measure = measure;
          best_parameters = grouprot.coordinate_system_parameters();
          if (verbosity > 1)
            std::cout << "new best " << best_measure << grouprot.coordinate_system_parameters() << std::endl;
        }
        if (measure < .1)
          break;
      }
      group.coordinate_system_parameters() = best_parameters;
      if (best_measure > 10)
        return false;
    }
    //  group.coordinate_system_parameters()[3]=-1.3796;
    //  group.coordinate_system_parameters()[4]=1.63646;
    //  group.coordinate_system_parameters()[5]= -2.13267;
    //    sm.reset_neighbours();
    //    std::cout << "chosen initial axes\n" << group.coordinate_system().axes() << std::endl;
    //    std::cout << "best scanned measure " << sm() << std::endl;
    //  for (int i=0; i<24; i++)
    //    std::cout << "group element "<<i<<", measure="<<sm(i,1,2)<<std::endl;
    //    std::cout << "Atomic coordinates in current local frame\n" << std::endl;
    //    for (const auto& atom : molecule.m_atoms)
    //      std::cout << atom.name << ": " << group.coordinate_system().to_local(atom.position).transpose() <<
    //      std::endl;
  } else { // not spherical_top()
    const std::string& s = group.name();
    std::smatch m;
    if (std::regex_match(s, m, std::regex{"Oh|O|Td|Th|T|Ih|I"}))
      return false;
  }
  //  if (sm() > threshold*1000) return false;
  SymmetryMeasure sm(molecule, group);
  sm.refine_frame(verbosity);
  double d = sm();
  if (verbosity >= 0)
    std::cout << "end of test_group() for " << group.name() << ", measure=" << d << " " << (d < threshold) << " "
              << threshold << std::endl;
  return d < threshold;
}
Eigen::Matrix3d find_axis_frame(const Molecule& molecule, const Group& group) {
  auto toprot = group.highest_rotation(true, 0);
  //  std::cout << "enter find_axis_frame, group " << group.name() << std::endl;
  //  std::cout << toprot << "\norder? " << toprot.order() << std::endl;
  if (toprot.order() < 3)
    return Eigen::Matrix3d::Identity();
  auto axis = molecule.findaxis(toprot.order());
  if (axis.norm() == 0)
    return Eigen::Matrix3d::Identity();
  //    throw std::runtime_error("failure to find molecular rotation axis of order " + std::to_string(toprot.order()));
//  std::cout << "found axis " << axis.transpose() << std::endl;
//  std::cout << "coordinate system rotation:\n" << group.coordinate_system().axes() << std::endl;
  //  auto toprotaxis =  group.coordinate_system().axes()*toprot.axis();
  auto toprotaxis = toprot.axis();
  //  std::cout << "toprot.axis "<< toprotaxis.transpose()<<std::endl;
  Eigen::Matrix3d R;
  R.col(2) = axis.cross(toprotaxis);
  //  std::cout << "r3 " << R.col(2).transpose() << std::endl;
  R.col(0) = -R.col(2).cross(axis);
  R.col(1) = R.col(2).cross(R.col(0));
  for (int i = 0; i < 3; i++)
    R.col(i) /= R.col(i).norm();
  //      std::cout << "rotation frame\n" << R << std::endl;
  //    std::cout << "determinant of rotation frame\n"<<R.determinant()<<std::endl;
  //    std::cout << "R(dag)*R\n" << R.transpose() * R <<std::endl;
  //      std::cout << "toprot.axis() "<<toprot.axis()<<" "<<axis.dot(toprot.axis()) <<" "<< axis.norm() <<" "<<
  //      toprot.axis().norm()<<std::endl;
  auto angle =
      std::acos(std::max(double(-1), std::min(double(1), axis.dot(toprotaxis) / (axis.norm() * toprotaxis.norm()))));
  //    std::cout << "angle=" << angle <<" rad = "<<angle*180/std::acos(double(-1))<< std::endl;

  if (std::abs(angle) < 1e-12)
    return Eigen::Matrix3d::Identity();
  Eigen::Matrix3d rot;
  rot << cos(angle), sin(angle), 0, -sin(angle), cos(angle), 0, 0, 0, 1;
  //    std::cout << "rotation in rotation frame\n"<<rot<<std::endl;
  //    std::cout << "determinant of rotation in rotation frame\n"<<rot.determinant()<<std::endl;
  auto new_axes = R * rot * R.inverse();
  //    std::cout << "determinant of new axes\n" << new_axes.determinant() << std::endl;
  //  std::cout << "new axes\n" << new_axes << std::endl;
  return new_axes;
}

static CoordinateSystem s_default_coordinate_system;

Group discover_group(const Molecule& molecule, double threshold, int verbosity) {
  return discover_group(molecule, s_default_coordinate_system, threshold, verbosity);
}

Group discover_group(const Molecule& molecule, CoordinateSystem& coordinate_system, double threshold, int verbosity) {
  auto p = molpro::Profiler::single()->push("SymmetryMeasure::discover_group(" + molecule.m_title + ")");
  using vec = CoordinateSystem::vec;
  const vec xaxis{1, 0, 0};
  const vec yaxis{0, 1, 0};
  const vec zaxis{0, 0, 1};
  const bool generators_only = true;
  constexpr size_t maximum_axis_order = 10;
  Group result;
  // special?
  for (const auto& n : std::vector<std::string>{"R3", "Dinfh", "Cinfv", "Oh", "O", "Td", "Ih", "I"})
    if (test_group(molecule, Group(coordinate_system, n, generators_only), threshold, verbosity))
      return Group(coordinate_system, n);

  // axis?
  for (int axis_order = maximum_axis_order; axis_order > 1; axis_order--) {
    auto o = std::to_string(axis_order);
    Group c2x(coordinate_system);
    c2x.name() = "pseudo-C2x";
    c2x.add(Rotation(coordinate_system, {0, 0, 1}, axis_order));
    double angle_break_symmetry = 0.01;
    c2x.add(Rotation(coordinate_system, {std::cos(angle_break_symmetry), std::sin(angle_break_symmetry), 0}, 2));
    Group sigma_h(coordinate_system);
    sigma_h.name() = "pseudo-sigma_h";
    sigma_h.add(Rotation(coordinate_system, {0, 0, 1}, axis_order));
    sigma_h.add(Reflection(coordinate_system, {0, 0, 1}));
    //    std::cout << "explore "
    //              << "C" + o << std::endl;
    //    std::cout << "Atomic coordinates in current local frame\n" << std::endl;
    //    for (const auto& atom : molecule.m_atoms)
    //      std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
    if (test_group(molecule, Group(coordinate_system, "C" + o, generators_only), threshold, verbosity)) {
      //      std::cout << "c2x\n" << c2x << std::endl;
      //      std::cout << "explore "
      //                << "c2x" << std::endl;
      //      std::cout << "Atomic coordinates in current local frame\n" << std::endl;
      //      for (const auto& atom : molecule.m_atoms)
      //        std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
      if (test_group(molecule, c2x, threshold, verbosity)) {
        //        std::cout << "test sigma_h: " << test_group(molecule, sigma_h, threshold,verbosity) << std::endl;
        //        if (test_group(molecule, sigma_h, threshold,verbosity)) {
        auto coordinate_system_save = coordinate_system;
        if (test_group(molecule, Group(coordinate_system, "D" + o + "h", generators_only), threshold, verbosity)) {
          return Group(coordinate_system, "D" + o + "h");
        } else {
          coordinate_system = coordinate_system_save;
          for (const auto& n : std::vector<std::string>{"D" + o + "d", "D" + o}) {
            if (verbosity > 1) {
              std::cout << "explore " << n << std::endl;
              std::cout << "Atomic coordinates in current local frame\n" << std::endl;
              for (const auto& atom : molecule.m_atoms)
                std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
            }
            if (test_group(molecule, Group(coordinate_system, n, generators_only), threshold, verbosity))
              return Group(coordinate_system, n);
          }
        }
      } else {
        auto o2 = std::to_string(axis_order * 2);
        for (const auto& n : std::vector<std::string>{"C" + o + "h", "C" + o + "v", "S" + o2, "C" + o}) {
          //       std::cout << "explore "<<n<<std::endl;
          if (n != "S2" and test_group(molecule, Group(coordinate_system, n, generators_only), threshold, verbosity))
            return Group(coordinate_system, n);
        }
      }
    }
  }
  // no axis found
  for (const auto& n : std::vector<std::string>{"Cs", "Ci", "C1"})
    if (test_group(molecule, Group(coordinate_system, n, generators_only), threshold, verbosity))
      return Group(coordinate_system, n);
  throw std::logic_error("unexpected failure to find point group");
}
Molecule molecule_localised(const CoordinateSystem& coordinate_system, const Molecule& source) {
  Molecule result(source);
  for (auto& atom : result.m_atoms) {
    atom.position = coordinate_system.to_local(atom.position);
  }
  return result;
}

} // namespace molpro::point_charge_symmetry
