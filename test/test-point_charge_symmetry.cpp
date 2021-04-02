#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iostream>
#include <molpro/Profiler.h>
#include <molpro/point_charge_symmetry/Group.h>
#include <molpro/point_charge_symmetry/Molecule.h>
#include <molpro/point_charge_symmetry/Operator.h>
#include <molpro/point_charge_symmetry/SymmetryMeasure.h>
#include <numeric>

using std::cout;
using namespace molpro::point_charge_symmetry;
using vec = Operator::vec;
using mat = Operator::mat;

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
  //    std::cout << "reference:\n" << reference << std::endl;
  for (int logstep = -7; logstep < 1; logstep++) {
    std::vector<mat> tested;
    const auto step = std::pow(double(10), logstep);
    for (int displacements = 1; displacements < 3; displacements++) {
      auto axes_gradient = coords.axes_gradient(displacements, step);
      tested.push_back(mat::Zero());
      for (int i = 0; i < 3; i++)
        tested.back() += displacement[i] * axes_gradient[i] / displacement.norm();
      //      std::cout << "tested:\n" << tested.back() << std::endl;
      //      std::cout << "reference-tested:\n" << reference - tested.back() << std::endl;
      const auto tolerance = std::max(1e-8, std::pow(step, displacements * 2));
      //      std::cout << "step=" << step << " , displacements=" << displacements
      //                << ", reference-tested=" << (reference - tested.back()).norm() << ", tolerance=" << tolerance
      //                << std::endl;
      EXPECT_LT((reference - tested.back()).norm(), tolerance);
    }
    const auto tolerance = std::max(1e-9, std::pow(step, 2));
    EXPECT_LT((tested[1] - tested[0]).norm(), tolerance);
  }
}

TEST(point_charge_symmetry, Molecule) {
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
  Molecule water("h2o-nosym.xyz");
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
  int i = 0;
  for (const auto &op : group) {

    std::cout << "Operator symmetry measure: " << op->name() << " " << sm(i) << std::endl;
    i++;
  }
  std::cout << sm << std::endl;
  std::cout << group.name() << " symmetry measure: " << sm() << std::endl;
  //  std::cout << "CoordinateSystem data";
  //  for (int i = 0; i < 6; i++)
  //    std::cout << " " << group.coordinate_system().data()[i];
  //  std::cout << std::endl;
  for (int functional_form = 0; functional_form < 2; functional_form++) {

    auto g = sm.coordinate_system_gradient(-1, functional_form);
    //    std::cout << "functional_form="<<functional_form<<std::endl;
    //    std::cout << "coordinate system gradient:";
    //    std::for_each(g.begin(), g.end(), [](const auto &val) { std::cout << " " << val; });
    //    std::cout << std::endl;

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
    //    std::cout << "numerical coordinate system gradient:";
    //    std::for_each(numerical_gradient.begin(), numerical_gradient.end(),
    //                  [](const auto &val) { std::cout << " " << val; });
    //    std::cout << std::endl;
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
  orders["C1"] = 1;
  orders["Ci"] = 2;
  orders["Cs"] = 2;
  orders["C2"] = 2;
  orders["C2v"] = 4;
  orders["S2"] = 2;
  orders["C2h"] = 4;
  orders["D2"] = 4;
  orders["D4"] = 8;
  orders["D2d"] = 8;
  orders["D4d"] = 16;
  orders["D2h"] = 8;
  orders["C3v"] = 6;
  orders["C3h"] = 6;
  orders["D3h"] = 12;
  orders["D3d"] = 12;
  orders["S4"] = 4;
  orders["S12"] = 12;
  for (const auto &n : orders) {
    auto g = molpro::point_charge_symmetry::group_factory(n.first);
    //    std::cout << g << std::endl;
    EXPECT_EQ(g.end() - g.begin(), n.second) << "Wrong order for group " << g;
  }
}

TEST(point_charge_symmetry, discover_group) {
  std::shared_ptr<molpro::Profiler> prof = molpro::Profiler::single("Discover groups");
  std::map<std::string, std::string> expected_groups;
  expected_groups["n2"] = "Dinfh";
  expected_groups["h2o"] = "C2v";
  expected_groups["h2o-nosym"] = "C2v";
  expected_groups["ferrocene"] = "D5d";
  expected_groups["benzene"] = "D6h";
  expected_groups["allene"] = "D2d";
  expected_groups["ch4"] = "Td";
  expected_groups["methane"] = "Td";
  expected_groups["p4"] = "Td";
  //  expected_groups["hexamethylbenzene"]="D3d";
  //  expected_groups["buckminsterfullerene"]="Ih";
  //  expected_groups["sulfur-hexafluoride"]="Oh";
  expected_groups["cyclohexane"] = "D3d";
  expected_groups["s8"] = "D4h";
  for (const auto &n : expected_groups) {
    Molecule molecule(n.first + ".xyz");
    auto group = molpro::point_charge_symmetry::discover_group(molecule, 1e-2);
    EXPECT_EQ(group.name(), n.second) << n.first << ": " << group.name();
    std::cout << n.first << ": " << group.name() << std::endl;
  }
  std::cout << *prof << std::endl;
}
