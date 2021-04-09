#include "SymmetryMeasure.h"
#include <molpro/Profiler.h>
#include <molpro/linalg/itsolv/OptimizeBFGS.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <regex>
#include <sstream>
#include <unsupported/Eigen/MatrixFunctions>
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

double SymmetryMeasure::operator()(int operator_index, int functional_form, int verbosity) const {
  auto p = molpro::Profiler::single()->push("SymmetryMeasure()");
  constexpr bool use_neighbour_list = true;
  double result = 0;
  auto start = operator_index < 0 ? m_group.begin() : m_group.begin() + operator_index;
  auto end = operator_index < 0 ? m_group.end() : m_group.begin() + operator_index + 1;
  if (verbosity > -1)
    std::cout << "!! SymmetryMeasure() " << std::endl;
  for (auto op = start; op < end; op++) {
    if (verbosity > 0) {
      std::cout << "Operator " << (*op)->name() << std::endl;
      std::cout << "Frame parameters";
      for (int i = 0; i < 6; i++)
        std::cout << " " << (*op)->coordinate_system().data()[i];
      std::cout << std::endl;
    }
    for (int a = 0; a < m_molecule.m_atoms.size(); a++) {
      int ai = m_neighbours[op - m_group.begin()][a];
      const auto& neighbour_position = m_molecule.m_atoms[ai].position;
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

void SymmetryMeasure::adopt_inertial_axes() {
  auto& coordinate_system = m_group.coordinate_system();
  auto& parameters = m_group.coordinate_system_parameters();
  //    std::cout << "adopt_inertial_axes initial, group="<<m_group.name()<<"\n" << coordinate_system.axes() <<
  //    std::endl;
  { // initialise to inertial axes in some orientation or other
    CoordinateSystem temporary_coordinate_system(m_molecule.centre_of_charge(), m_molecule.inertial_axes());
    std::copy(temporary_coordinate_system.m_parameters.begin(), temporary_coordinate_system.m_parameters.end(),
              parameters.begin());
  }
  //  std::cout << "adopt_inertial_axes first inertial\n" << coordinate_system.axes() << std::endl;
  //  std::cout << "centre of charge " << m_molecule.centre_of_charge() << std::endl;
  //  std::cout <<"Group: "<<m_group<<std::endl;
  int best_axis = 0;
  double best_axis_sm = 1e50;
  bool symmetric_top = false;
  for (int principal_axis = 0; principal_axis < 6; principal_axis++) {
    reset_neighbours();
    //        std::cout << "try axes\n" << coordinate_system.axes() << std::endl;
    //    auto sm = SymmetryMeasure(m_molecule, m_group);
    auto local_inertia_tensor =
        (coordinate_system.axes().transpose() * m_molecule.inertia_tensor() * coordinate_system.axes()).eval();
    //    std::cout << "local_inertia_tensor\n" << local_inertia_tensor << std::endl;
    double measure;
    m_inertia_principal_values = local_inertia_tensor.diagonal().eval();
    if (std::abs(local_inertia_tensor(0, 0) - local_inertia_tensor(1, 1)) <
        1e-5 * std::abs(local_inertia_tensor(0, 0) - local_inertia_tensor(2, 2))) {
      //      std::cout << local_inertia_tensor(0, 0) << " " << local_inertia_tensor(1, 1) << std::endl;
      //      std::cout << "Symmetric top with axes\n" << coordinate_system.axes() << std::endl;
      measure = 0;
      symmetric_top = true;
    } else
      measure = (*this)();
    //        std::cout << "Atomic coordinates in local frame\n" << std::endl;
    //    for (const auto &atom : m_molecule.m_atoms)
    //      std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
    CoordinateSystem::vec out_of_plane{0, 0, 0};
    for (const auto& atom : m_molecule.m_atoms) {
      //      auto coords= coordinate_system.to_local(atom.position);
      out_of_plane += coordinate_system.to_local(atom.position).cwiseAbs();
      //      std::cout << coords.transpose()<<std::endl;
      //      std::cout << coordinate_system.to_local(atom.position).cwiseAbs().transpose()<<std::endl;
      //      std::cout << out_of_plane.transpose()<<std::endl;
      //      std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
    }
    // TODO implement all recommended tie breakers, such as planar C2v molecules in yz plane
    auto tie_breaker = out_of_plane.dot(CoordinateSystem::vec{3, 2, 1}) * 1e-6;
    measure += tie_breaker;
    //    std::cout << "principal_axis=" << principal_axis << " sm=" << measure << std::endl;
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
  //      std::cout << "axes before symmetric top adjustment\n" << coordinate_system.axes() << std::endl;
  if (false and symmetric_top) { // add random rotation to step off possible maximum
                                 //    coordinate_system.m_parameters[3] += 0*std::acos(double(-1))/4;
    auto oldaxes = coordinate_system.axes();
    auto sq2 = std::sqrt(double(0.5));
    CoordinateSystem::mat transform;
    transform << sq2, sq2, 0, -sq2, sq2, 0, 0, 0, 1;
    coordinate_system.from_axes((oldaxes * transform).eval());
  }
  reset_neighbours();
  //  std::cout << "chosen initial coordinate system: " << best_axis << "\n" << coordinate_system << std::endl;
  //  std::cout << "Atomic coordinates in local frame\n" << std::endl;
  //  for (const auto& atom : m_molecule.m_atoms)
  //    std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
}

int SymmetryMeasure::optimise_frame() {
  constexpr bool optimize_origin = false;
  auto& parameters = m_group.coordinate_system_parameters();
  const double centre_of_charge_penalty = 0e0;
  const int centre_of_charge_penalty_power = 4;
  const int verbosity = -1;
  using Rvector = CoordinateSystem::parameters_t;
  //  std::cout << "optimise_frame" << std::endl;
  for (int c = 0; c < 0; c++) { // TODO remove testing only
    parameters = {1, 1, 1, 1, 1, 1};
    double step = 1e-3;
    auto value0 = (*this)();
    auto grad0 = coordinate_system_gradient();
    parameters[c] += step;
    auto valuep = (*this)();
    std::cout << "plus displacement parameters=";
    for (int i = 0; i < 6; i++)
      std::cout << " " << parameters[i];
    std::cout << ", value=" << valuep << std::endl;
    parameters[c] -= 2 * step;
    auto valuem = (*this)();
    std::cout << "minus displacement parameters=";
    for (int i = 0; i < 6; i++)
      std::cout << " " << parameters[i];
    std::cout << ", value=" << valuem << std::endl;
    parameters[c] += step;
    std::cout << "analytic=" << grad0[c] << ", numerical=" << (valuep - valuem) / (2 * step) << std::endl;
  }
  auto solver = molpro::linalg::itsolv::create_Optimize<Rvector, Rvector, Rvector>(
      "BFGS", "max_size_qspace=3,convergence_threshold=1e-6");
  int nwork = 1;
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
  for (int iter = 0; iter < 200; iter++) {
    //    std::cout << coordinate_system << std::endl;
    reset_neighbours();
    auto value = (*this)(-1, 1);
    auto grad = coordinate_system_gradient(-1, 1);
    if (false) { // test gradient
      CoordinateSystem::parameters_t vpp, vp, vm, vmm, testgrad;
      double step = 2e-3;
      double diff = 0;
      for (int i = 0; i < 6; i++) {
        parameters[i] -= 2 * step;
        vmm[i] = (*this)(-1, 1);
        parameters[i] += step;
        vm[i] = (*this)(-1, 1);
        parameters[i] += 2 * step;
        vp[i] = (*this)(-1, 1);
        parameters[i] += step;
        vpp[i] = (*this)(-1, 1);
        parameters[i] -= 2 * step;
        testgrad[i] = (vmm[i] - 8 * vm[i] + 8 * vp[i] - vpp[i]) / (12 * step);
        diff += std::pow(testgrad[i] - grad[i], 2);
      }
      //      std::cout << "calc gradient " << grad << std::endl;
      //      std::cout << "test gradient " << testgrad << std::endl;
      diff = std::sqrt(diff);
      //      std::cout << "diff gradient " << diff << std::endl;
      if (diff > 1e-7)
        throw std::runtime_error("differentiation");
    }
    if (not optimize_origin)
      std::fill(grad.begin(), grad.begin() + 3, 0);
    auto centre_of_charge_displacement = (m_group.coordinate_system().origin() - m_molecule.centre_of_charge()).eval();
    if (verbosity > 1) {
      std::cout << "grad " << grad << std::endl;
      std::cout << "coordinate_system.origin() " << m_group.coordinate_system().origin().transpose() << std::endl;
      std::cout << "m_molecule.centre_of_charge() " << m_molecule.centre_of_charge().transpose() << std::endl;
      std::cout << "centre_of_charge_displacement " << centre_of_charge_displacement.transpose() << std::endl;
    }
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
    //    std::cout << "before add_value " << nwork;
    //    for (int i = 0; i < 6; i++)
    //      std::cout << " " << parameters[i];
    //    std::cout << std::endl;
    //    std::cout << " before add_value " << nwork;
    //    for (int i = 0; i < 6; i++)
    //      std::cout << " " << grad[i];
    //    std::cout << std::endl;
    //    std::cout << "value=" << value << std::endl;
    //    std::cout << "Current parameters " << parameters << std::endl;
    auto current_parameters = parameters;
    //    grad[0] = 0;
    //    grad[1] = 0;
    //    grad[2] = 0;
    if (solver->add_value(parameters, value, grad)) {
      if (verbosity > 1) {
        std::cout << "after add_value " << nwork;
        for (int i = 0; i < 6; i++)
          std::cout << " " << grad[i];
        std::cout << std::endl;
        cout << "precondition" << std::endl;
      }
      for (int i = 0; i < 6; i++)
        grad[i] /= 100;
    } else if (false) {
      std::cout << "LINE SEARCH WAS SPECIFIED" << std::endl;
      auto new_parameters = parameters;
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
          parameters[i] = (1 - x) * current_parameters[i] + x * new_parameters[i];
        int component = 5;
        auto fi = (*this)(-1, 1);
        parameters[component] += 1e-2;
        auto fid = (*this)(-1, 1);
        parameters[component] -= 1e-2;
        auto gradi = coordinate_system_gradient(-1, 1);
        std::cout << "component " << component << " numerical gradient=" << (fi - fid) / 1e-2
                  << ", analytical gradient=" << gradi[component] << std::endl;
        double grad_proj = 0;
        for (int i = 0; i < 6; i++)
          grad_proj += gradi[i] * (new_parameters[i] - current_parameters[i]);
        std::cout << "x=" << x << " value=" << (*this)(-1, 1) << ", projected gradient " << grad_proj << std::endl;
      }
      parameters = new_parameters;
    }
    nwork = solver->end_iteration(parameters, grad);
    if (verbosity > 1) {
      std::cout << "after end_iteration " << nwork;
      for (int i = 0; i < 6; i++)
        std::cout << " " << parameters[i];
      std::cout << std::endl;
      std::cout << m_group.coordinate_system().axes() << std::endl;
      std::cout << "Atomic coordinates in current local frame\n" << std::endl;
      for (const auto& atom : m_molecule.m_atoms)
        std::cout << atom.name << ": " << m_group.coordinate_system().to_local(atom.position).transpose() << std::endl;
    }
    if (verbosity > 0)
      solver->report();
    if (nwork <= 0)
      return iter;
  }
  return -1;
}

Molecule SymmetryMeasure::refine() const {
  auto result = m_molecule;
  return result;
}

static inline bool test_group(const Molecule& molecule, const Group& group, double threshold = 1e-6,
                              int verbosity = -1) {
  if (verbosity >= 0)
    std::cout << "test_group " << group.name() << std::endl;
  if (verbosity > 0)
    std::cout << "test_group " << group << std::endl;
  auto p = molpro::Profiler::single()->push("SymmetryMeasure::test_group(" + group.name() + ")");
  SymmetryMeasure sm(molecule, group);
  sm.adopt_inertial_axes();
  if (verbosity > 0) {
    std::cout << "initial measure " << sm() << std::endl;
    std::cout << group.coordinate_system().axes() << std::endl;
    std::cout << "Atomic coordinates in current local frame\n" << std::endl;
  }
  // scan
  if (sm.spherical_top()) {
    const int nscan = 5;
    double best_measure = 1e50;
    double pi = std::acos(double(-1));
    CoordinateSystem::parameters_t best_parameters;
    for (int zscan = 0; zscan < nscan; zscan++)
      for (int yscan = 0; yscan < nscan; yscan++)
        for (int xscan = 0; xscan < nscan; xscan++) {
          group.coordinate_system_parameters()[3] = (2 * zscan + 1 - nscan) * pi / nscan;
          group.coordinate_system_parameters()[4] = (2 * yscan + 1 - nscan) * pi / nscan;
          group.coordinate_system_parameters()[5] = (2 * xscan + 1 - nscan) * pi / nscan;
          sm.reset_neighbours();
          auto measure = sm();
          //        std::cout << "try "<<measure<<" ( current best="<<best_measure<<")
          //        "<<group.coordinate_system_parameters()<<std::endl;
          if (measure < best_measure) {
            best_measure = measure;
            best_parameters = group.coordinate_system_parameters();
            //            std::cout << "new best " << best_measure << group.coordinate_system_parameters() << std::endl;
          }
        }
    group.coordinate_system_parameters() = best_parameters;
    //  group.coordinate_system_parameters()[3]=-1.3796;
    //  group.coordinate_system_parameters()[4]=1.63646;
    //  group.coordinate_system_parameters()[5]= -2.13267;
    sm.reset_neighbours();
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
  sm.optimise_frame();
  double d = sm();
  if (verbosity >= 0)
    std::cout << "end of test_group() for " << group.name() << ", measure=" << d << " " << (d < threshold) << " "
              << threshold << std::endl;
  return d < threshold;
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
  const bool generators_only = false;
  constexpr size_t maximum_axis_order = 6;
  Group result;
  using molpro::point_charge_symmetry::group_factory;
  // special?
  //  for (const auto& n : std::vector<std::string>{"Dinfh", "Cinfv", "Td", "Oh", "Ih"})
  //  for (const auto& n : std::vector<std::string>{"Oh", "Td"})
  for (const auto& n : std::vector<std::string>{"Dinfh", "Cinfv", "Td"})
    if (test_group(molecule, group_factory(coordinate_system, n, generators_only), threshold, verbosity))
      return group_factory(coordinate_system, n);

  // axis?
  for (int axis_order = maximum_axis_order; axis_order > 1; axis_order--) {
    auto o = std::to_string(axis_order);
    Group c2x(coordinate_system, "pseudo-C2x");
    c2x.add(Rotation(coordinate_system, {0, 0, 1}, axis_order));
    c2x.add(Rotation(coordinate_system, {std::cos(double(.001)), std::sin(double(.001)), 0}, 2));
    Group sigma_h(coordinate_system, "pseudo-sigma_h");
    sigma_h.add(Rotation(coordinate_system, {0, 0, 1}, axis_order));
    sigma_h.add(Reflection(coordinate_system, {0, 0, 1}));
    //    std::cout << "explore "
    //              << "C" + o << std::endl;
    //    std::cout << "Atomic coordinates in current local frame\n" << std::endl;
    //    for (const auto& atom : molecule.m_atoms)
    //      std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
    if (test_group(molecule, group_factory(coordinate_system, "C" + o, generators_only), threshold, verbosity)) {
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
        if (test_group(molecule, group_factory(coordinate_system, "D" + o + "h", generators_only), threshold,
                       verbosity)) {
          return group_factory(coordinate_system, "D" + o + "h");
        } else {
          coordinate_system = coordinate_system_save;
          for (const auto& n : std::vector<std::string>{"D" + o + "d", "D" + o}) {
            if (verbosity > 1) {
              std::cout << "explore " << n << std::endl;
              std::cout << "Atomic coordinates in current local frame\n" << std::endl;
              for (const auto& atom : molecule.m_atoms)
                std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
            }
            if (test_group(molecule, group_factory(coordinate_system, n, generators_only), threshold, verbosity))
              return group_factory(coordinate_system, n);
          }
        }
      } else {
        for (const auto& n : std::vector<std::string>{"C" + o + "h", "C" + o + "v", "S" + o, "C" + o}) {
          //       std::cout << "explore "<<n<<std::endl;
          if (n != "S2" and
              test_group(molecule, group_factory(coordinate_system, n, generators_only), threshold, verbosity))
            return group_factory(coordinate_system, n);
        }
      }
    }
  }
  // no axis found
  for (const auto& n : std::vector<std::string>{"Cs", "Ci", "C1"})
    if (test_group(molecule, group_factory(coordinate_system, n, generators_only), threshold, verbosity))
      return group_factory(coordinate_system, n);
  throw std::logic_error("unexpected failure to find point group");
}

Group group_factory(std::string name, bool generators_only) {
  return group_factory(s_default_coordinate_system, name, generators_only);
}
Group group_factory(CoordinateSystem& coordinate_system, std::string name, bool generators_only) {
  using vec = CoordinateSystem::vec;
  const vec xaxis{1, 0, 0};
  const vec yaxis{0, 1, 0};
  const vec zaxis{0, 0, 1};
  const auto all = not generators_only;
  Group g(coordinate_system, name);
  if (all)
    g.add(Identity());
  std::smatch m;
  if (name == "Ih")
    g.add(Rotation(zaxis, 7)); // TODO implement
  if (name == "Oh") {
    for (int axis = 0; axis < 3; axis++) {
      for (int count = 0; count < 3; count += 2) {
        g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, true, count));
        if (all)
          g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, false, count));
      }
      g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 2, true, 1));
      g.add(Reflection(Eigen::Matrix3d::Identity().col(axis)));
    }
    for (int corner = 0; corner < (all ? 4 : 2); corner++) {
      auto sq2 = std::sqrt(1 / double(2));
      g.add(Rotation(vec{std::cos((2 * corner + 1) * acos(double(-1)) / 4),
                         std::sin((2 * corner + 1) * acos(double(-1)) / 4), sq2},
                     3));
    }
  }
  if (name == "Td") {
    for (int axis = 0; axis < 3; axis++) {
      if (all)
        g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 2, true, 1));
      if (all)
        g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, false, 1));
      if (all)
        g.add(Rotation(Eigen::Matrix3d::Identity().col(axis), 4, false, 3));
      if (all)
        g.add(Reflection(Eigen::Matrix3d::Identity().col((axis + 1) % 3) +
                         Eigen::Matrix3d::Identity().col((axis + 2) % 3)));
      g.add(Reflection(Eigen::Matrix3d::Identity().col((axis + 1) % 3) -
                       Eigen::Matrix3d::Identity().col((axis + 2) % 3)));
    }
    for (int corner = 0; corner < (all ? 4 : 2); corner++) {
      auto sq2 = std::sqrt(1 / double(2));
      for (int count = 1; count < 3; count++)
        g.add(Rotation(vec{std::cos((2 * corner + 1) * acos(double(-1)) / 4),
                           std::sin((2 * corner + 1) * acos(double(-1)) / 4), sq2},
                       3, true, count));
    }
  }
  if (name == "Cinfv" or name == "Dinfh") { // representative only
    auto g = group_factory(coordinate_system, name.replace(1, 3, "11"), generators_only);
    g.name() = name;
    return g;
  }

  if (std::regex_match(name, m, std::regex{"[CD][0-9]*[02468]h"}) or
      std::regex_match(name, m, std::regex{"D[0-9]*[13579]d"}) or
      std::regex_match(name, m, std::regex{"Ih|Th|Oh|Dinfh|Ci"}))
    g.add(Inversion());

  if (std::regex_match(name, m, std::regex{"Dinfh|Ih|Cs|[CD][1-9][0-9]*h"}))
    g.add(Reflection(zaxis));

  if (std::regex_match(name, m, std::regex{"([C])([1-9][0-9]*)([v])"}) or
      std::regex_match(name, m, std::regex{"([D])([1-9][0-9]*)([hd])"})) {
    auto order = std::stoi(m.str(2));
    for (double angle = 0; angle < std::acos(double(-1)) - 1e-10; angle += std::acos(double(-1)) / order)
      g.add(Reflection({std::cos(angle), std::sin(angle), 0}));
  }

  if (std::regex_match(name, m, std::regex{"([D])([1-9][0-9]*)([^v])?"})) {
    auto order = std::stoi(m.str(2));
    //    for (double angle = 0; angle < std::acos(double(-1)) - 1e-10; angle += std::acos(double(-1)) / order)
    for (int count = 0; count < (all ? order : 2); count++) {
      double angle = (std::regex_match(name,m,std::regex{"D[0-9]*[02468]d"}) ? (count + 0.5) : count) * std::acos(double(-1)) / order;
      g.add(Rotation({std::cos(angle), std::sin(angle), 0}, 2));
    }
  }

  if (std::regex_match(name, m, std::regex{"([CD])([1-9][0-9]*)([hvd]*)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 1; count < order; count++)
      g.add(Rotation(zaxis, order, true, count));
  }
  if (std::regex_match(name, m, std::regex{"([CD])([0-9]*[02468])(h)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 0; count < order / 2; count++)
      if (count * 2 + 1 != order / 2)
        g.add(Rotation(zaxis, order, false, count * 2 + 1));
    if (order > 2)
      for (int count = 0; count < ((order % 4) ? order / 2 : order / 4); count++)
        if (count * 2 + 1 != ((order % 4) ? order / 2 : order / 4))
          g.add(Rotation(zaxis, order / 2, false, count * 2 + 1));
  }
  if (std::regex_match(name, m, std::regex{"([CD])([0-9]*[13579])(h)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 0; count < order; count++)
      if (count != (order - 1) / 2)
        g.add(Rotation(zaxis, order, false, count * 2 + 1));
  }
  if (std::regex_match(name, m, std::regex{"([D])([1-9][0-9]*)(d)"})) {
    auto order = std::stoi(m.str(2));
    for (int count = 0; count < order; count++)
      if (2 * count != (order - 1))
        g.add(Rotation(zaxis, 2 * order, false, count * 2 + 1));
  }
  if (std::regex_match(
          name, m,
          std::regex{
              "([S])([0-9]*[02468])"})) { // could be done more prettily with explicit proper rotations and inversion
    auto order = std::stoi(m.str(2));
    for (int count = 1; count < order; count++)
      g.add(Rotation(zaxis, order, false, count));
  }
  return g;
}

Molecule molecule_localised(const CoordinateSystem& coordinate_system, const Molecule& source) {
  Molecule result (source);
  for (auto& atom : result.m_atoms) {
    atom.position = coordinate_system.to_local(atom.position);
  }
  return result;
}

} // namespace molpro::point_charge_symmetry
