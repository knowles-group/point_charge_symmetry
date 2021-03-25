#include "CoordinateSystem.h"
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

namespace molpro::point_charge_symmetry {

CoordinateSystem::CoordinateSystem(const vec& origin, const mat& axes) {
  if (std::abs(axes.determinant() - 1) > 1e-14)
    throw std::runtime_error("Axes must be proper rotation");
  auto generator = axes.log();
  axis_generator()(0) = generator(1, 0);
  axis_generator()(1) = generator(2, 0);
  axis_generator()(2) = generator(2, 1);
  this->origin()(0) = origin(0);
  this->origin()(1) = origin(1);
  this->origin()(2) = origin(2);
}

const CoordinateSystem::mat CoordinateSystem::axes() const {
  Eigen::Matrix3d generator;
  generator << 0, axis_generator()(0), axis_generator()(1), -axis_generator()(0), 0, axis_generator()(2),
      -axis_generator()(1), -axis_generator()(2), 0;
  //  std::cout << "Generator\n"<<generator<<std::endl;
  return generator.exp();
}

const std::array<CoordinateSystem::mat, 3> CoordinateSystem::axes_gradient(int displacements, double step) const {
  std::array<CoordinateSystem::mat, 3> result;
  for (int i = 0; i < 3; i++) {
    std::vector<CoordinateSystem> points;
    for (int displacement = -displacements; displacement <= displacements; displacement++) {
      points.push_back(*this);
      points.back().axis_generator()[i] += step * displacement;
    }
    if (displacements == 1)
      result[i] = (points[2].axes() - points[0].axes()) / (2 * step);
    else if (displacements == 2)
      result[i] = (-points[4].axes() + 8 * points[3].axes() + points[0].axes() - 8 * points[1].axes()) / (12 * step);
    else
      throw std::logic_error("Incorrect differentiation order");
  }
  return result;
}
} // namespace molpro::point_charge_symmetry
