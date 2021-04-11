#include <array>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iostream>
#include <molpro/Profiler.h>
#include <molpro/point_charge_symmetry/Group.h>
#include <molpro/point_charge_symmetry/Molecule.h>
#include <molpro/point_charge_symmetry/Operator.h>
#include <molpro/point_charge_symmetry/SymmetryMeasure.h>
#include <numeric>
#include <vector>

using std::cout;
using namespace molpro::point_charge_symmetry;
using vec = Operator::vec;
using mat = Operator::mat;

TEST(point_charge_symmetry, Euler) {
  CoordinateSystem cs;
  cs.m_rotation_parameter_type = RotationParameterType::Euler;
  std::vector<std::array<double, 3>> parameters;
  parameters.push_back({1, 0, 0});
  parameters.push_back({0, 1, 0});
  parameters.push_back({0, 0, 1});
  parameters.push_back({0, 1e-13, 1});
  parameters.push_back({0, 1e-12, 1});
  parameters.push_back({0, 1e-11, 1});
  parameters.push_back({0, 1e-10, 1});
  parameters.push_back({0, 1e-9, 1});
  parameters.push_back({0, 1e-8, 1});
  parameters.push_back({0, 1e-7, 1});
  parameters.push_back({0, 1e-6, 1});
  parameters.push_back({1, 1, 0});
  parameters.push_back({0.1, 0.2, 0.3});
  for (const auto &p : parameters) {
    //    cout << "test "<<p[0]<<", "<<p[1]<<", "<<p[2]<<std::endl;
    for (int i = 0; i < 3; i++)
      cs.m_parameters[3 + i] = p[i];
    //    cout << cs << std::endl;
    auto u = cs.axes();
    cs.from_axes(u);
    auto pnew = std::array<double, 3>{cs.m_parameters[3], cs.m_parameters[4], cs.m_parameters[5]};
    if (std::abs(p[1]) > 1e-12)
      EXPECT_THAT(pnew, ::testing::Pointwise(::testing::DoubleNear(1e-10), p));
    //    cout << cs << std::endl;
    auto unew = cs.axes();
    EXPECT_THAT(std::vector<double>(&unew(0, 0), &unew(0, 0) + 9),
                ::testing::Pointwise(::testing::DoubleNear(1e-10), std::vector<double>(&u(0, 0), &u(0, 0) + 9)));
  }
}

void test_operation(const vec &initial, const Operator &op, const vec &expected) {
  auto ttt = op(initial);
  EXPECT_THAT(
      std::vector<double>(ttt.data(), ttt.data() + 3),
      ::testing::Pointwise(::testing::DoubleNear(1e-13), std::vector<double>(expected.data(), expected.data() + 3)))
      << "original: " << initial.transpose() << "\nOperator: " << op;
}

TEST(point_charge_symmetry, local_operations) {
  test_operation({1, 1, 1}, Reflection({0, 0, 1}), {1, 1, -1});
  test_operation({1, 1, 1}, Reflection({0, 0, -1}), {1, 1, -1});
  test_operation({1, 1, 1}, Reflection({-1, -1, -1}), {-1, -1, -1});
  test_operation({1, 1, 1}, Reflection({1, -1, 0}), {1, 1, 1});
  test_operation({1, 1, 1}, Reflection({-1, 1, 0}), {1, 1, 1});
  test_operation({-1, -1, 1}, Reflection({-1, 1, 0}), {-1, -1, 1});
  test_operation({1, 1, 1}, Rotation({0, 0, 1}, 2), {-1, -1, 1});
  test_operation({1, 1, 1}, Rotation({0, 0, 1}, 4), {-1, 1, 1});
  test_operation({1, 1, 1}, Rotation({0, 0, 1}, 4, false), {-1, 1, -1});
  test_operation({1, 1, 1}, Inversion(), {-1, -1, -1});
}

TEST(point_charge_symmetry, translated_operations) {
  CoordinateSystem coords({10, 10, 10});
  test_operation({1, 1, 1}, Reflection(coords, {0, 0, 1}), {1, 1, 19});
  test_operation({1, 1, 1}, Reflection(coords, {0, 0, -1}), {1, 1, 19});
  test_operation({1, 1, 1}, Reflection(coords, {-1, -1, -1}), {19, 19, 19});
  test_operation({1, 1, 1}, Reflection(coords, {1, -1, 0}), {1, 1, 1});
  test_operation({1, 1, 1}, Reflection(coords, {-1, 1, 0}), {1, 1, 1});
  test_operation({-1, -1, 1}, Reflection(coords, {-1, 1, 0}), {-1, -1, 1});
}

TEST(point_charge_symmetry, rotated_operations) {
  mat axes;
  axes << 0, 1, 0, -1, 0, 0, 0, 0, 1;
  CoordinateSystem coords({0, 0, 0}, axes);
  test_operation({1, 1, 1}, Reflection(coords, {0, 0, 1}), {1, 1, -1});
  test_operation({1, 1, 1}, Reflection(coords, {0, 0, -1}), {1, 1, -1});
  test_operation({1, 1, 1}, Reflection(coords, {-1, -1, -1}), {1 / 3., 5 / 3., 1 / 3.});
  test_operation({1, 1, 1}, Reflection(coords, {1, -1, 0}), {-1, -1, 1});
  test_operation({1, 1, 1}, Reflection(coords, {-1, 1, 0}), {-1, -1, 1});
  test_operation({-1, -1, 1}, Reflection(coords, {-1, 1, 0}), {1, 1, 1});
}

TEST(point_charge_symmetry, Group) {
  mat axes;
  axes << 0, 1, 0, -1, 0, 0, 0, 0, 1;
  CoordinateSystem coords({0, 0, 0});
  Group c2v(coords);
  c2v.add(Identity());
  c2v.add(Rotation({0, 0, 1}, 2));
  c2v.add(Reflection({1, 0, 0}));
  c2v.add(Reflection({0, 1, 0}));
}

TEST(point_charge_symmetry, axes_gradient) {
  std::shared_ptr<molpro::Profiler> prof = molpro::Profiler::single("axes_gradient");
  mat axes;
  //  axes << 0, 1, 0, -1, 0, 0, 0, 0, 1;
  axes << 1 / std::sqrt(3), 1 / std::sqrt(3), 1 / std::sqrt(3), 2 / std::sqrt(6), -1 / std::sqrt(6), -1 / std::sqrt(6),
      0, 1 / std::sqrt(2), -1 / std::sqrt(2);
  CoordinateSystem coords({0, 0, 0}, axes);
  //  std::cout << "coords.axes():\n" << coords.axes() << std::endl;
  vec displacement{2e-6, 2e-6, 2e-6};
  auto coords_minus = coords;
  coords_minus.axis_generator() -= displacement;
  auto coords_plus = coords;
  coords_plus.axis_generator() += displacement;
  auto reference = ((coords_plus.axes() - coords_minus.axes()) / (2 * displacement.norm())).eval();
  //  std::cout << "coords.axes():\n" << coords.axes() << std::endl;
  //  std::cout << "coords_plus.axes():\n" << coords_plus.axes() << std::endl;
  //  std::cout << "coords_minus.axes():\n" << coords_minus.axes() << std::endl;
  //      std::cout << "reference:\n" << reference << std::endl;
  for (int logstep = -7; logstep < 1; logstep++) {
    std::vector<mat> tested;
    const auto step = std::pow(double(10), logstep);
    for (int displacements = 0; displacements < 3; displacements++) {
      auto axes_gradient = coords.axes_gradient(displacements, step);
      tested.emplace_back(mat::Zero());
      for (int i = 0; i < 3; i++)
        tested.back() += displacement[i] * axes_gradient[i] / displacement.norm();
      //            std::cout << "tested:\n" << tested.back() << std::endl;
      //            std::cout << "reference-tested:\n" << reference - tested.back() << std::endl;
      const auto tolerance = std::max(1e-8, 2 * std::pow(step, displacements * 2));
      EXPECT_LT((reference - tested.back()).norm(), tolerance)
          << "step=" << step << " , displacements=" << displacements
          << ", reference-tested=" << (reference - tested.back()).norm() << ", tolerance=" << tolerance << std::endl;
    }
    const auto tolerance = std::max(2e-9, 2 * std::pow(step, 2));
    EXPECT_LT((tested[1] - tested[0]).norm(), tolerance);
  }
}

TEST(point_charge_symmetry, Molecule) {
  std::shared_ptr<molpro::Profiler> prof = molpro::Profiler::single("Molecule");
  Molecule water("h2o.xyz");
  std::cout << water << std::endl;
  Group c2v("C2v");
  c2v.add(Identity());
  c2v.add(Rotation({0, 0, 1}, 2));
  c2v.add(Reflection({1, 0, 0}));
  c2v.add(Reflection({0, 1, 0}));
  auto sm = SymmetryMeasure(water, c2v);
  std::cout << sm << std::endl;
  int i = 0;
  for (const auto &op : c2v) {

    std::cout << "Operator symmetry measure: " << op->name() << " " << sm(i) << std::endl;
    i++;
  }
  std::cout << c2v.name() << " symmetry measure: " << sm() << std::endl;
  //  std::cout << "CoordinateSystem data";
  //  for (int i = 0; i < 6; i++)
  //    std::cout << " " << c2v.coordinate_system().data()[i];
  //  std::cout << std::endl;
}

TEST(point_charge_symmetry, SymmetryMeasure_gradient) {
  std::shared_ptr<molpro::Profiler> prof = molpro::Profiler::single("SymmetryMeasure_gradient");
  Molecule water("h2o-nosym.xyz");
  //  Molecule water("h2o.xyz");
  //  Molecule water("Ferrocene.xyz");
  std::cout << water << std::endl;
  //  std::cout << "centre of charge: " << water.centre_of_charge().transpose() << std::endl;
  std::cout << "inertial axes\n" << water.inertial_axes() << std::endl;
  auto axes = water.inertial_axes();
  //  CoordinateSystem coordinate_system(water.centre_of_charge(), axes);
  CoordinateSystem coordinate_system;
  //  Group group(coordinate_system, "special");
  //  group.add(Rotation({0, 0, 1}, 5,true,1));
  //  group.add(Rotation({0, 0, 1}, 3,true,1));
  //  group.add(Identity());
  //  group.add(Rotation({0, 0, 1}, 2));
  //  group.add(Reflection({1, 0, 0}));
  //  group.add(Reflection({0, 1, 0}));
  auto group = group_factory(coordinate_system, "C3v");
  //  int best_axis = 0;
  //  double best_axis_sm = 1e50;
  //  for (int principal_axis = 0; principal_axis < 6; principal_axis++) {
  //        std::cout << "try axes\n" << coordinate_system.axes() << std::endl;
  //    auto sm = SymmetryMeasure(water, group);
  ////    std::cout << "Atomic coordinates in local frame\n" << std::endl;
  ////    for (const auto &atom : water.m_atoms)
  ////      std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
  //    std::cout << "principal_axis=" << principal_axis << " sm=" << sm() << std::endl;
  //    if (sm() < best_axis_sm) {
  //      best_axis_sm=sm();
  //      best_axis = principal_axis;
  //    }
  //    //    std::cout << "axes before cycle\n" << coordinate_system.axes() << std::endl;
  //    coordinate_system.cycle_axes();
  //  }
  //  for (int principal_axis = 0; principal_axis < best_axis; principal_axis++) {
  //    coordinate_system.cycle_axes();
  //  }
  //  std::cout << "chosen initial coordinate system: " << best_axis << "\n" << coordinate_system << std::endl;
  auto sm = SymmetryMeasure(water, group);
  sm.adopt_inertial_axes();
  std::cout << "Atomic coordinates in local frame\n" << std::endl;
  for (const auto &atom : water.m_atoms)
    std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
  std::cout << sm << std::endl;
  std::cout << group.name() << " symmetry measure: " << sm() << std::endl;
  //  std::cout << "CoordinateSystem data";
  //  for (int i = 0; i < 6; i++)
  //    std::cout << " " << group.coordinate_system().data()[i];
  //  std::cout << std::endl;
  for (int functional_form = 0; functional_form < 2; functional_form++) {

    auto g = sm.coordinate_system_gradient(-1, functional_form);
    std::cout << "functional_form=" << functional_form << std::endl;
    std::cout << "coordinate system gradient:";
    std::for_each(g.begin(), g.end(), [](const auto &val) { std::cout << " " << val; });
    std::cout << std::endl;

    std::array<double, 6> numerical_gradient{0, 0, 0, 0, 0, 0};
    for (int i = 0; i < 6; i++) {
      //      std::cout << "Displace " << i << std::endl;
      double *parameters = coordinate_system.data();
      //      std::fill(parameters, parameters + 6, 0);
      double step = 1e-6;
      parameters[i] += step;
      //    std::cout << "CoordinateSystem data";
      //    for (int i = 0; i < 6; i++)
      //      std::cout << " " << group.coordinate_system().data()[i];
      //    std::cout << std::endl;
      auto smp = sm(-1, functional_form);
      //          std::cout << smp<< std::endl;
      parameters[i] -= 2 * step;
      auto smm = sm(-1, functional_form);
      parameters[i] += step;
      numerical_gradient[i] = (smp - smm) / (2 * step);
      //      std::cout << "displaced " << i << " smm=" << smm << ", smp=" << smp << std::endl;
    }
    std::cout << "numerical coordinate system gradient:";
    std::for_each(numerical_gradient.begin(), numerical_gradient.end(),
                  [](const auto &val) { std::cout << " " << val; });
    std::cout << std::endl;
    EXPECT_THAT(numerical_gradient, ::testing::Pointwise(::testing::DoubleNear(1e-8), g))
        << "for functional_form=" << functional_form;
  }

  EXPECT_GE(sm.optimise_frame(), 0);
  EXPECT_GE(sm.optimise_frame(), 0);
  std::cout << group.name() << " symmetry measure: " << sm() << std::endl;
  std::cout << coordinate_system << std::endl;
  std::cout << "Atomic coordinates in local frame\n" << std::endl;
  for (const auto &atom : water.m_atoms)
    std::cout << atom.name << ": " << coordinate_system.to_local(atom.position).transpose() << std::endl;
}

TEST(point_charge_symmetry, group_factory) {
  //  std::cout<<"D5h generators\n" << molpro::point_charge_symmetry::group_factory("D5h",true)<<std::endl;
  std::map<std::string, int> orders;
  orders["Dinfh"] = 44;
  orders["Cinfv"] = 22;
  orders["C1"] = 1;
  orders["Ci"] = 2;
  orders["Cs"] = 2;
  orders["Td"] = 24;
  //  orders["Th"] = 24;
  //  orders["T"] = 12;
  //  orders["Oh"] = 48;
  //  orders["O"] = 24;
  //  orders["Ih"] = 120;
  //  orders["I"] = 60;
  for (int i = 2; i < 12; i++) {
    orders[std::string{"S"} + std::to_string(2 * i)] = i * 2;
    orders[std::string{"D"} + std::to_string(i) + "h"] = i * 4;
    orders[std::string{"D"} + std::to_string(i) + "d"] = i * 4;
    orders[std::string{"D"} + std::to_string(i)] = i * 2;
    orders[std::string{"C"} + std::to_string(i) + "h"] = i * 2;
    orders[std::string{"C"} + std::to_string(i) + "v"] = i * 2;
    orders[std::string{"C"} + std::to_string(i)] = i;
  }
  for (const auto &n : orders) {
    auto g = molpro::point_charge_symmetry::group_factory(n.first);
    //    std::cout << g << std::endl;
    EXPECT_EQ(g.end() - g.begin(), n.second) << "Wrong order for group " << g;
  }
}

TEST(point_charge_symmetry, discover_group) {
  std::shared_ptr<molpro::Profiler> prof = molpro::Profiler::single("Discover groups");
  prof->set_max_depth(1);
  std::map<std::string, std::string> expected_groups;
  //  expected_groups["n2"] = "Dinfh";
  expected_groups["h2o"] = "C2v";
  expected_groups["h2o-nosym"] = "C2v";
  expected_groups["ferrocene"] = "D5d";
  expected_groups["benzene"] = "D6h";
  expected_groups["allene"] = "D2d";
  expected_groups["ch4"] = "Td";
  expected_groups["ethane"] = "D3d";
  expected_groups["methane"] = "Td";
  expected_groups["p4"] = "Td";
  expected_groups["adamantane"] = "Td";
  //  expected_groups["hexamethylbenzene"]="D3d";
  //  expected_groups["buckminsterfullerene"]="Ih";
  //  expected_groups["sulfur-hexafluoride"]="Oh";
  expected_groups["cyclohexane"] = "D3d";
  expected_groups["s8"] = "D8h";
  expected_groups["phloroglucinol"] = "C3h";
  for (const auto &n : expected_groups) {
    Molecule molecule(n.first + ".xyz");
    auto group = molpro::point_charge_symmetry::discover_group(molecule, 1e-2, -1);
    EXPECT_EQ(group.name(), n.second) << n.first << ": " << group.name();
    std::cout << n.first << ": " << group.name() << ", measure=" << SymmetryMeasure(molecule, group)() << std::endl;
  }
  std::cout << *prof << std::endl;
}
TEST(point_charge_symmetry, allene45) {
  std::shared_ptr<molpro::Profiler> prof = molpro::Profiler::single("allene45");
  Molecule allene("allene45.xyz");
  //  std::cout << allene<<std::endl;
  CoordinateSystem::mat axes;
  axes << 0, 0, 1, 1, 0, 0, 0, 1, 0;
  CoordinateSystem cs(vec::Zero(), axes);
  //  std::cout << cs<<std::endl;
  const Group group = group_factory(cs, "D2d");
  //  std::cout << group << std::endl;
  SymmetryMeasure sm(allene, group);
  EXPECT_LE(sm(-1, 0, -1), 1e-14);
  auto grad = sm.coordinate_system_gradient();
  EXPECT_LE(std::inner_product(grad.begin(), grad.end(), grad.begin(), 0), 1e-14); //<< "Gradient: " << grad;
}

template <typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &v) {
  for (const auto &e : v)
    s << " " << e;
  return s;
}

TEST(point_charge_symmetry, atom_gradient) {
  std::shared_ptr<molpro::Profiler> prof = molpro::Profiler::single("atom_gradient");
  Molecule molecule("hexamethylbenzene.xyz");
  CoordinateSystem cs;
  //  auto group = discover_group(molecule, cs);
  auto group = group_factory(cs, "D3h");
  SymmetryMeasure sm(molecule, group);
  //  std::cout << molecule << std::endl;
  //  std::cout << sm() << std::endl;
  auto grad = sm.atom_gradient();
  //  cout << "analytic gradient " << grad << std::endl;
  auto numgrad = grad;
  numgrad.assign(numgrad.size(), 0);
  int i = 0;
  double step = 1e-3;
  auto sm0 = sm();
  for (auto &atom : molecule.m_atoms) {
    for (int j = 0; j < 3; j++) {
      atom.position(j) -= 2 * step;
      auto smmm = sm();
      atom.position(j) += step;
      auto smm = sm();
      atom.position(j) += 2 * step;
      auto smp = sm();
      atom.position(j) += step;
      auto smpp = sm();
      atom.position(j) -= 2 * step;
      numgrad[i] = (smmm - 8 * smm + 8 * smp - smpp) / (12 * step);
      i++;
    }
  }
  //  cout << "numerical gradient " << numgrad << std::endl;
  EXPECT_THAT(grad, ::testing::Pointwise(::testing::DoubleNear(1e-6), numgrad));
}
