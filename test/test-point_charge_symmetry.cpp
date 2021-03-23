#include <gtest/gtest.h>
#include <iostream>
#include <molpro/point_charge_symmetry/Operator.h>
#include <gmock/gmock.h>
#include <molpro/point_charge_symmetry/Group.h>
#include <molpro/point_charge_symmetry/Molecule.h>

using std::cout;
using namespace molpro::point_charge_symmetry;
using vec = Operator::vec;
using mat = Operator::mat;

void test_operation(const vec &initial, const Operator &op, const vec &expected) {
  auto ttt = op(initial);
  EXPECT_THAT(std::vector<double>(ttt.data(),ttt.data()+3),
  ::testing::Pointwise(::testing::DoubleNear(1e-13), std::vector<double>(expected.data(),expected.data()+3)))<< "original: "<<initial.transpose() << "\nOperator: "<<op;
}
TEST(point_charge_symmetry, local_operations) {
  test_operation({1,1,1},ReflectionPlane({0,0,1}),{1,1,-1});
  test_operation({1,1,1},ReflectionPlane({0,0,-1}),{1,1,-1});
  test_operation({1,1,1},ReflectionPlane({-1,-1,-1}),{-1,-1,-1});
  test_operation({1,1,1},ReflectionPlane({1,-1,0}),{1,1,1});
  test_operation({1,1,1},ReflectionPlane({-1,1,0}),{1,1,1});
  test_operation({-1,-1,1},ReflectionPlane({-1,1,0}),{-1,-1,1});
  test_operation({1,1,1},Axis({0,0,1},2),{-1,-1,1});
  test_operation({1,1,1},Axis({0,0,1},4),{-1,1,1});
  test_operation({1,1,1},Axis({0,0,1},4,false),{-1,1,-1});
  test_operation({1,1,1},Inversion(),{-1,-1,-1});
}

TEST(point_charge_symmetry, translated_operations) {
  CoordinateSystem coords({10,10,10});
  test_operation({1,1,1},ReflectionPlane(coords,{0,0,1}),{1,1,19});
  test_operation({1,1,1},ReflectionPlane(coords,{0,0,-1}),{1,1,19});
  test_operation({1,1,1},ReflectionPlane(coords,{-1,-1,-1}),{19,19,19});
  test_operation({1,1,1},ReflectionPlane(coords,{1,-1,0}),{1,1,1});
  test_operation({1,1,1},ReflectionPlane(coords,{-1,1,0}),{1,1,1});
  test_operation({-1,-1,1},ReflectionPlane(coords,{-1,1,0}),{-1,-1,1});
}

TEST(point_charge_symmetry, rotated_operations) {
  mat axes; axes << 0,1,0,-1,0,0,0,0,1;
  CoordinateSystem coords({0,0,0},axes);
  test_operation({1,1,1},ReflectionPlane(coords,{0,0,1}),{1,1,-1});
  test_operation({1,1,1},ReflectionPlane(coords,{0,0,-1}),{1,1,-1});
  test_operation({1,1,1},ReflectionPlane(coords,{-1,-1,-1}),{5/3.,1/3.,1/3.});
  test_operation({1,1,1},ReflectionPlane(coords,{1,-1,0}),{-1,-1,1});
  test_operation({1,1,1},ReflectionPlane(coords,{-1,1,0}),{-1,-1,1});
  test_operation({-1,-1,1},ReflectionPlane(coords,{-1,1,0}),{1,1,1});
}

TEST(point_charge_symmetry, Group) {
  mat axes; axes << 0,1,0,-1,0,0,0,0,1;
  CoordinateSystem coords({0,0,0});
  Group c2v(coords);
  c2v.add(Identity());
  c2v.add(Axis({0,0,1},2));
  c2v.add(ReflectionPlane({1,0,0}));
  c2v.add(ReflectionPlane({0,1,0}));
}

TEST(point_charge_symmetry, Molecule) {
  Molecule water("h2o.xyz");
  std::cout << water << std::endl;
}