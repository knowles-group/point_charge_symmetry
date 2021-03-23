#include <iostream>
#include <molpro/point_charge_symmetry/SymmetryOperator.h>

int main() {
  using namespace molpro::point_charge_symmetry;
  ReflectionPlane thing({1,0,0});
  std::cout << "Hello, World!" << std::endl;
  return 0;
}
