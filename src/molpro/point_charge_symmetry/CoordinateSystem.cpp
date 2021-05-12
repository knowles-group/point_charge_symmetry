#include "CoordinateSystem.h"
#include <iostream>
#include <sstream>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <molpro/Profiler.h>

namespace molpro::point_charge_symmetry {

CoordinateSystem::CoordinateSystem(const vec& origin, const mat& axes) {
  //  std::cout << "CoordinateSystem constructor, axes=\n" << axes << std::endl;
  this->origin()(0) = origin(0);
  this->origin()(1) = origin(1);
  this->origin()(2) = origin(2);
  from_axes(axes);
}

void CoordinateSystem::from_axes(const mat& axes) const {
  if (std::abs(axes.determinant() - 1) > 1e-14)
    throw std::runtime_error("Axes must be proper rotation");
  switch (m_rotation_parameter_type) {
  case RotationParameterType::Log: {
    auto generator = axes.log();
    m_parameters[3 + 2] = generator(1, 0);
    m_parameters[3 + 1] = generator(2, 0);
    m_parameters[3 + 0] = generator(2, 1);
  } break;
  case RotationParameterType::Euler: {
    const auto& u = axes;
    auto& alpha = m_parameters[3 + 0];
    auto& beta = m_parameters[3 + 1];
    auto& gamma = m_parameters[3 + 2];
    beta = std::acos(u(2, 2));
    if (beta < 1e-6)
      beta = std::asin(std::sqrt(std::pow(u(0, 2), 2) + std::pow(u(1, 2), 2)));
    if (beta < 1e-12) {
      alpha = std::atan2(u(1, 0), u(0, 0));
      gamma = 0;
    } else {
      alpha = std::atan2(u(1, 2), u(0, 2));
      gamma = std::atan2(u(2, 1), -u(2, 0));
      auto pi = std::acos(double(-1));
      if (gamma >= pi)
        gamma -= pi;
    }
  } break;
  case RotationParameterType::Quaternion:
    throw std::logic_error("unimplemented");
    break;
  }
}

const CoordinateSystem::mat CoordinateSystem::axes() const {
  switch (m_rotation_parameter_type) {
  case RotationParameterType::Log: {
    Eigen::Matrix3d generator;
    generator << 0, -axis_generator()(2), -axis_generator()(1), axis_generator()(2), 0, -axis_generator()(0),
        axis_generator()(1), axis_generator()(0), 0;
    //    std::cout << "Generator\n"<<generator<<std::endl;
    //  std::cout << "Generator.exp()\n"<<generator.exp()<<std::endl;
    return generator.exp();
  } break;
  case RotationParameterType::Euler: {
    auto c1 = std::cos(axis_generator()(0));
    auto s1 = std::sin(axis_generator()(0));
    auto c2 = std::cos(axis_generator()(1));
    auto s2 = std::sin(axis_generator()(1));
    auto c3 = std::cos(axis_generator()(2));
    auto s3 = std::sin(axis_generator()(2));
    Eigen::Matrix3d result;
    result << c1 * c2 * c3 - s1 * s3, -c3 * s1 - c1 * c2 * s3, c1 * s2, c2 * c3 * s1 + c1 * s3, c1 * c3 - c2 * s1 * s3,
        s1 * s2, -c3 * s2, s2 * s3, c2;
    return result;
  } break;
  case RotationParameterType::Quaternion:
    throw std::logic_error("unimplemented");
    break;
  }
  throw std::logic_error("unexpected rotation parameter type");
}

const std::array<CoordinateSystem::mat, 3> CoordinateSystem::axes_gradient(int displacements, double step) const {
  std::array<CoordinateSystem::mat, 3> result;
//  auto p = molpro::Profiler::single()->push("CoordinateSystem::axes_gradient()");
  if (displacements == 0) { // analytic
    if (m_rotation_parameter_type == RotationParameterType::Euler) {
      auto c1 = std::cos(axis_generator()(0));
      auto s1 = std::sin(axis_generator()(0));
      auto c2 = std::cos(axis_generator()(1));
      auto s2 = std::sin(axis_generator()(1));
      auto c3 = std::cos(axis_generator()(2));
      auto s3 = std::sin(axis_generator()(2));
      result[0] << -c2 * c3 * s1 - c1 * s3, -c1 * c3 + c2 * s1 * s3, -s1 * s2, c1 * c2 * c3 - s1 * s3,
          -c3 * s1 - c1 * c2 * s3, c1 * s2, 0, 0, 0;
      result[1] << -c1 * c3 * s2, c1 * s2 * s3, c1 * c2, -c3 * s1 * s2, s1 * s2 * s3, c2 * s1, -c2 * c3, c2 * s3, -s2;
      result[2] << -c3 * s1 - c1 * c2 * s3, -c1 * c2 * c3 + s1 * s3, 0, c1 * c3 - c2 * s1 * s3, -c2 * c3 * s1 - c1 * s3,
          0, s2 * s3, c3 * s2, 0;
    } else
      throw std::logic_error("Unimplemented");
  } else {
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
  }
  return result;
}

CoordinateSystem::vec CoordinateSystem::to_local(const vec& source) const {
  return (axes().transpose() * (source - origin())).eval();
}

CoordinateSystem::vec CoordinateSystem::to_global(const vec& source) const {
  return (origin() + axes() * source).eval();
}

void CoordinateSystem::cycle_axes() const {
  mat new_axes;
  mat axes = this->axes();
  if (m_axis_permutation_rot90_next) {
    new_axes.col(0) = axes.col(1);
    new_axes.col(1) = -axes.col(0);
    new_axes.col(2) = axes.col(2);
  } else {
    new_axes.col(0) = axes.col(2);
    new_axes.col(1) = -axes.col(1);
    new_axes.col(2) = axes.col(0);
  }
  for (int i = 0; i < 2; i++)
    for (int j = i + 1; j < 3; j++)
      if (new_axes(i, i) < 0 and new_axes(j, j) < 0) {
        new_axes.col(i) = -new_axes.col(i);
        new_axes.col(j) = -new_axes.col(j);
      }
  from_axes(new_axes);
  m_axis_permutation_rot90_next = not m_axis_permutation_rot90_next;
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
