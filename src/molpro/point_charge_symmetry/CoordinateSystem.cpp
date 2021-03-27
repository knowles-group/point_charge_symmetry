#include "CoordinateSystem.h"
#include <iostream>
#include <sstream>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

namespace molpro::point_charge_symmetry {

CoordinateSystem::CoordinateSystem(const vec& origin, const mat& axes) {
  std::cout << "CoordinateSystem constructor, axes=\n"<<axes<<std::endl;
  if (std::abs(axes.determinant() - 1) > 1e-14)
    throw std::runtime_error("Axes must be proper rotation");
  auto generator = axes.log();
  axis_generator()(2) = generator(1, 0);
  axis_generator()(1) = generator(2, 0);
  axis_generator()(0) = generator(2, 1);
  this->origin()(0) = origin(0);
  this->origin()(1) = origin(1);
  this->origin()(2) = origin(2);
}

const CoordinateSystem::mat CoordinateSystem::axes() const {
  Eigen::Matrix3d generator;
  generator << 0, -axis_generator()(2), -axis_generator()(1), axis_generator()(2), 0, -axis_generator()(0),
      axis_generator()(1), axis_generator()(0), 0;
//    std::cout << "Generator\n"<<generator<<std::endl;
//  std::cout << "Generator.exp()\n"<<generator.exp()<<std::endl;
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

CoordinateSystem::vec CoordinateSystem::to_local(const vec& source) {
  return (axes().transpose() * (source - origin())).eval();
}

CoordinateSystem::vec CoordinateSystem::to_global(const vec& source) { return (origin() + axes() * source).eval(); }

void CoordinateSystem::cycle_axes() {
  mat new_axes;
  mat axes = this->axes();
//  std::cout << "cycle_axes() start\n"<<axes<<std::endl;
  new_axes.col(0) = axes.col(2);
  new_axes.col(1) = axes.col(0);
  new_axes.col(2) = axes.col(1);
//  std::cout << "cycle_axes() new\n"<<new_axes<<std::endl;
  auto generator = new_axes.log();
//  std::cout << "generator\n"<<generator<<std::endl;
  axis_generator()(2) = generator(1, 0);
  axis_generator()(1) = generator(2, 0);
  axis_generator()(0) = generator(2, 1);
}

void CoordinateSystem::rot90(int axis) {
  mat new_axes;
  mat axes = this->axes();
  new_axes.col(axis) = axes.col(axis);
  new_axes.col((axis + 1) % 3) = axes.col((axis + 2) % 3);
  new_axes.col((axis + 2) % 3) = -axes.col((axis + 1) % 3);
  auto generator = new_axes.log();
  axis_generator()(2) = generator(1, 0);
  axis_generator()(1) = generator(2, 0);
  axis_generator()(0) = generator(2, 1);
}

std::string CoordinateSystem::str() const {
  std::stringstream ss;
  ss << "CoordinateSystem origin =";
  for (int i = 0; i < 3; i++)
    ss << " " << m_parameters[i];
  ss << " axis generators =";
  for (int i = 0; i < 3; i++)
    ss << " " << m_parameters[i + 3];
  ss << "\nAxes:\n" << axes();
  return ss.str();
}
} // namespace molpro::point_charge_symmetry
