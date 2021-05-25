#include "Operator.h"
//#include <iostream>
#include <memory>
//#include <molpro/Profiler.h>
#include <sstream>
#include <vector>

namespace molpro::point_charge_symmetry {

static CoordinateSystem s_default_coordinate_system = CoordinateSystem();

Operator::Operator() : Operator(s_default_coordinate_system) {}

Operator::vec Operator::operator()(vec v) const {
  return m_coordinate_system.origin() +
         m_coordinate_system.axes() *
             operator_local((m_coordinate_system.axes().transpose() * (v - m_coordinate_system.origin()))).eval();
}
std::array<Operator::vec, 6> Operator::operator_gradient(vec v, int numerical, double step) const {
  //  auto p = molpro::Profiler::single()->push("Operator::operator_gradient()");
  std::array<Operator::vec, 6> result;
  if (numerical > 0) {
    for (int i = 0; i < 6; i++) {
      std::vector<std::unique_ptr<Operator>> points;
      std::vector<CoordinateSystem> coordinate_systems;
      std::vector<vec> transformed_points;
      for (int displacement = -numerical; displacement <= numerical; displacement++) {
        coordinate_systems.emplace_back(m_coordinate_system);
        coordinate_systems.back().data()[i] += step * displacement;
        points.emplace_back(this->clone(coordinate_systems.back()));
        transformed_points.push_back((*points.back())(v).eval());
        //        std::cout << "transformed_points.back() "<<transformed_points.back().transpose()<<std::endl;
      }
      if (numerical == 1)
        result[i] = (transformed_points[2] - transformed_points[0]) / (2 * step);
      else if (numerical == 2)
        result[i] =
            (transformed_points[0] - 8 * transformed_points[1] + 8 * transformed_points[3] - transformed_points[4]) /
            (12 * step);
      else
        throw std::logic_error("Incorrect differentiation order");
    }
  } else { // analytic differentation
    auto u = m_coordinate_system.axes();
    auto du = m_coordinate_system.axes_gradient();
    auto vshift = v - m_coordinate_system.origin();
    for (int i = 0; i < 3; i++) {
      result[i] = mat::Identity().col(i) - u * operator_local(u.transpose() * mat::Identity().col(i));
      result[3 + i] = u * operator_local(du[i].transpose() * vshift) + du[i] * operator_local(u.transpose() * vshift);
    }
  }
  //  std::cout << "result ";
  //  for (int i = 0; i < 6; i++)
  //    std::cout << " {" << result[i].transpose()<<"}";
  //  std::cout << std::endl;
  return result;
}

Reflection::Reflection(vec normal) : Reflection(s_default_coordinate_system, std::move(normal)) {}
Reflection::Reflection(const CoordinateSystem& coordinate_system, vec normal)
    : Operator(coordinate_system), m_normal(normal.normalized()) {
  m_name = "sigma";
  if (std::abs(m_normal(0)) > 0.99)
    m_name += "_yz";
  else if (std::abs(m_normal(1)) > 0.99)
    m_name += "_xz";
  else if (std::abs(m_normal(2)) > 0.99)
    m_name += "_h";
  else if (std::abs(m_normal(1)) < 0.01)
    m_name += "_v";
}
Operator::vec Reflection::operator_local(vec v) const {
  v -= 2 * m_normal.dot(v) * m_normal;
  return v;
}
Rotation::Rotation(vec axis, int order, bool proper, int count)
    : Rotation(s_default_coordinate_system, std::move(axis), std::move(order), std::move(proper), std::move(count)) {}
Rotation::Rotation(const CoordinateSystem& coordinate_system, vec axis, int order, bool proper, int count)
    : Operator(coordinate_system), m_axis(axis.normalized()), m_order(std::move(order)), m_proper(std::move(proper)),
      m_count(std::move(count)) {
  m_name = (m_proper ? "C" : "S") + std::to_string(m_order);
  if (m_axis(0) > 0.99)
    m_name += "x";
  if (m_axis(1) > 0.99)
    m_name += "y";
  if (m_axis(2) > 0.99)
    m_name += "z";
  if (m_count > 1)
    m_name += "^" + std::to_string(m_count);
}
Operator::vec Rotation::operator_local(vec v) const {
  auto aa = Eigen::AngleAxis<double>((double)2 * m_count * std::acos(double(-1)) / m_order, m_axis);
  v = aa * v;
  if (not m_proper)
    v -= 2 * m_axis.dot(v) * m_axis;
  return v;
}

Inversion::Inversion() : Operator(s_default_coordinate_system) {}
Inversion::Inversion(const CoordinateSystem& coordinate_system) : Operator(coordinate_system) { m_name = "i"; }
Operator::vec Inversion::operator_local(vec v) const { return -v; }

Identity::Identity() : Identity(s_default_coordinate_system) {}
Identity::Identity(const CoordinateSystem& coordinate_system) : Operator(coordinate_system) { m_name = "E"; }
Operator::vec Identity::operator_local(vec v) const { return v; }

std::string Operator::str(const std::string& title, bool coordinate_frame) const {
  std::stringstream result;
  result << m_name;
  if (!title.empty())
    result << " " << title;
  if (coordinate_frame) {
    result << ", origin: " << this->m_coordinate_system.origin().transpose();
    result << ", axes:\n" << this->m_coordinate_system.axes();
  }
  return result.str();
}

std::string Reflection::str(const std::string& title, bool coordinate_frame) const {
  std::stringstream result;
  result << "Reflection " << title;
  result << Operator::str(title, coordinate_frame);
  result << ", normal: " << this->m_normal.transpose();
  return result.str();
}

std::string Rotation::str(const std::string& title, bool coordinate_frame) const {
  std::stringstream result;
  result << (m_proper ? "R" : "Improper r") << "otation " << title;
  result << Operator::str(title, coordinate_frame);
  result << ", axis: " << this->m_axis.transpose();
  result << ", angle: " << this->m_count * double(360) / m_order;
  return result.str();
}

std::string Inversion::str(const std::string& title, bool coordinate_frame) const {
  std::stringstream result;
  result << "Inversion " << title;
  result << Operator::str(title, coordinate_frame);
  return result.str();
}
std::string Identity::str(const std::string& title, bool coordinate_frame) const {
  std::stringstream result;
  result << "Identity " << title;
  result << Operator::str(title, coordinate_frame);
  return result.str();
}
} // namespace molpro::point_charge_symmetry
