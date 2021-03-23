#include "CoordinateSystem.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

namespace molpro::point_charge_symmetry {

CoordinateSystem::CoordinateSystem(const vec &origin, const mat &axes) {
  if (std::abs(axes.determinant() - 1) > 1e-14) throw std::runtime_error("Axes must be proper rotation");
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
  generator << 0, axis_generator()(0), axis_generator()(1),
      -axis_generator()(0), 0, axis_generator()(2),
      -axis_generator()(1), -axis_generator()(2), 0;
  return generator.exp();
}
}
