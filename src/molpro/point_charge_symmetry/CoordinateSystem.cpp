#include "CoordinateSystem.h"
#include <iostream>
//#include <molpro/Profiler.h>
#include <cmath>
#include <sstream>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include "Euler.h"

namespace molpro::point_charge_symmetry {

CoordinateSystem::CoordinateSystem(const RotationParameterType rotation_parameter_type, const vec& origin,
                                   const mat& axes)
    : m_rotation_parameter_type(rotation_parameter_type) {
  //  std::cout << "CoordinateSystem constructor, axes=\n" << axes << std::endl;
  this->origin()(0) = origin(0);
  this->origin()(1) = origin(1);
  this->origin()(2) = origin(2);
  from_axes(axes);
}

CoordinateSystem& CoordinateSystem::operator=(const CoordinateSystem& source) {
  m_parameters = source.m_parameters;
  assert(m_rotation_parameter_type == source.m_rotation_parameter_type);
  m_axis_permutation_rot90_next = source.m_axis_permutation_rot90_next;
  return *this;
}

void CoordinateSystem::from_axes(const mat& axes) const {
  if (std::abs(axes.determinant() - 1) > 1e-14)
    throw std::runtime_error("Axes must be proper rotation");
  switch (m_rotation_parameter_type) {
  case RotationParameterType::Log: {
    auto generator = axes.log().eval();
    m_parameters[3 + 2] = generator(1, 0);
    m_parameters[3 + 1] = generator(2, 0);
    m_parameters[3 + 0] = generator(2, 1);
  } break;
  case RotationParameterType::Euler: {
    auto euler = euler_from_axes(axes);
    std::copy(euler.begin(), euler.end(), m_parameters.begin() + 3);
  } break;
  case RotationParameterType::TanEuler: {
    auto euler = euler_from_axes(axes);
    auto pi = std::acos(double(-1));
    m_parameters[3] = std::tan((euler[0] - pi) / (2 * m_rotation_parameter_scale));
    m_parameters[4] = std::tan(euler[1] / m_rotation_parameter_scale);
    m_parameters[5] = std::tan((euler[2] - pi) / (2 * m_rotation_parameter_scale));
    //        std::cout << "from axes\n" << axes << std::endl;
    //        std::cout << "Euler angles " << euler[0] << ", " << euler[1] << ", " << euler[2] << std::endl;
    //        std::cout << "parameters set to " << m_parameters[3] << ", " << m_parameters[4] << ", " << m_parameters[5]
    //                  << std::endl;
  } break;
  case RotationParameterType::Quaternion:
    throw std::logic_error("unimplemented");
    break;
  }
}
bool CoordinateSystem::operator==(const CoordinateSystem& other) const { return true; } // TODO implement
const CoordinateSystem::mat CoordinateSystem::axes() const {
  switch (m_rotation_parameter_type) {
  case RotationParameterType::Log: {
    Eigen::Matrix3d generator;
    generator << 0, -axis_generator()(2), -axis_generator()(1), axis_generator()(2), 0, -axis_generator()(0),
        axis_generator()(1), axis_generator()(0), 0;
    //    std::cout << "Generator\n"<<generator<<std::endl;
    //  std::cout << "Generator.exp()\n"<<generator.exp()<<std::endl;
    return generator.exp();
  }
  case RotationParameterType::Euler:
    return axes_from_euler(axis_generator());
  case RotationParameterType::TanEuler: {
    Eigen::Vector3d euler;
    auto pi = std::acos(double(-1));
    euler[0] = pi + 2 * m_rotation_parameter_scale * std::atan(axis_generator()(0));
    euler[1] = m_rotation_parameter_scale * std::atan(axis_generator()(1));
    euler[2] = pi + 2 * m_rotation_parameter_scale * std::atan(axis_generator()(2));
    //        std::cout << "parameters read " << m_parameters[3] << ", " << m_parameters[4] << ", " << m_parameters[5]
    //                  << std::endl;
    //        std::cout << "Euler angles " << euler[0] << ", " << euler[1] << ", " << euler[2] << std::endl;
    //        std::cout << "to axes\n" << axes_from_euler(euler) << std::endl;
    return axes_from_euler(euler);
  }
  case RotationParameterType::Quaternion:
    throw std::logic_error("unimplemented");
  }
  throw std::logic_error("unexpected rotation parameter type");
}
const std::array<CoordinateSystem::mat, 3> CoordinateSystem::axes_gradient(int displacements, double step) const {
  //  auto p = molpro::Profiler::single()->push("CoordinateSystem::axes_gradient()");
  if (displacements == 0) { // analytic
    switch (m_rotation_parameter_type) {
    case RotationParameterType::Log:
      throw std::logic_error("Unimplemented");
    case RotationParameterType::Euler:
      return axes_gradient_from_euler(axis_generator());
    case RotationParameterType::TanEuler: {
      std::array<CoordinateSystem::mat, 3> result;
      auto eulergrad = axes_gradient_from_euler(axis_generator());
      result[0] = 2 * m_rotation_parameter_scale * eulergrad[0] / (1 + std::pow(m_parameters[3], 2));
      result[1] = m_rotation_parameter_scale * eulergrad[1] / (1 + std::pow(m_parameters[4], 2));
      result[2] = 2 * m_rotation_parameter_scale * eulergrad[2] / (1 + std::pow(m_parameters[5], 2));
      return result;
    }
    case RotationParameterType::Quaternion:
      throw std::logic_error("Unimplemented");
    }
  } else {
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
  throw std::logic_error("unexpected");
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
std::array<std::array<double, 2>, 3> CoordinateSystem::rotation_generator_ranges() const {
  std::array<std::array<double, 2>, 3> result;
  const auto pi = std::acos(double(-1));
  switch (m_rotation_parameter_type) {
  case RotationParameterType::Log:
    result[0] = result[1] = result[2] = {0, 2 * pi};
    return result;
  case RotationParameterType::Euler:
    result[0] = result[2] = {0, 2 * pi};
    result[1] = {-pi / 2, pi / 2};
    return result;
  case RotationParameterType::TanEuler:;
  case RotationParameterType::Quaternion:;
  }
  throw std::logic_error("unimplemented");
}
} // namespace molpro::point_charge_symmetry
