#include "Operator.h"
#include <iostream>
#include <memory>
//#include <molpro/Profiler.h>
#include "Euler.h"
#include <sstream>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

namespace molpro::point_charge_symmetry {

static CoordinateSystem s_default_coordinate_system = CoordinateSystem();

Operator::Operator() : Operator(s_default_coordinate_system) {}

Operator::vec Operator::operator()(vec v) const {
  return m_coordinate_system.origin() +
         m_coordinate_system.axes() *
             operator_local((m_coordinate_system.axes().transpose() * (v - m_coordinate_system.origin()))).eval();
}
Operator::mat Operator::operator()() const {
  mat result;
  mat identity{mat::Identity()};
  for (Eigen::Index i = 0; i < 3; ++i)
    result.col(i) = (*this)(identity.col(i));
  return result;
}
Operator::vec Operator::operator_local(vec v) const { return m_local_representation * v; }
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

GenericOperator::GenericOperator() : GenericOperator(s_default_coordinate_system) {}
GenericOperator::GenericOperator(const CoordinateSystem& coordinate_system) : Operator(coordinate_system) {
  m_name = "generic";
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
  this->m_local_representation = mat::Identity() - m_normal * m_normal.transpose() * 2;
}
// Operator::vec Reflection::operator_local(vec v) const {
//   v -= 2 * m_normal.dot(v) * m_normal;
//   return v;
// }
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
  this->m_local_representation =
      Eigen::AngleAxis<double>((double)2 * m_count * std::acos(double(-1)) / m_order, m_axis);
  if (not m_proper)
    this->m_local_representation -= 2 * m_axis * m_axis.transpose();
}

Inversion::Inversion() : Inversion(s_default_coordinate_system) {}
Inversion::Inversion(const CoordinateSystem& coordinate_system) : Operator(coordinate_system) {
  m_name = "i";
  this->m_local_representation = -mat::Identity();
}

Identity::Identity() : Identity(s_default_coordinate_system) {}
Identity::Identity(const CoordinateSystem& coordinate_system) : Operator(coordinate_system) {
  m_name = "E";
  this->m_local_representation = mat::Identity();
}
GenericOperator::GenericOperator(const CoordinateSystem& coordinate_system, const GenericOperator& source)
    : Operator(coordinate_system) {
  this->m_name = source.m_name;
  this->m_local_representation = source.m_local_representation;
}
GenericOperator::GenericOperator(const Operator& A, const Operator& B) : Operator(A.coordinate_system()) {
  //  std::cout << "A: " << A << std::endl;
  //  std::cout << "B: " << B << std::endl;
  if (A.coordinate_system() != B.coordinate_system())
    throw std::runtime_error("Incompatible coordinate systems");
  m_name = A.name() + " * " + B.name();
  //  std::cout << "m_name " << m_name << std::endl;
  //  std::cout << "A():\n"<<A()<<std::endl;
  //  std::cout << "B():\n" << B() << std::endl;
  const mat BB = B().eval();
  const mat AA = A().eval();
  mat AB = (AA * BB).eval();
  for (Eigen::Index i = 0; i < 3; ++i)
    for (Eigen::Index j = 0; j < 3; ++j) {
      if (std::abs(AA(i, j)) > 1)
        throw std::runtime_error("unexpected values in A");
      if (std::abs(BB(i, j)) > 1)
        throw std::runtime_error("unexpected values in B");
      if (std::abs(AB(i, j)) > 1 + 1e-12)
        throw std::runtime_error("unexpected values in AB");
      AB(i, j) = std::min(double(1), std::max(AB(i, j), double(-1)));
      if (std::abs(AB(i, j)) < 1e-12)
        AB(i, j) = 0;
    }
  //  std::cout << "AB\n" << AB << std::endl;
  this->m_local_representation = AB.eval();
  //  m_local_representation = (m_local_representation * mat::Identity()).eval(); // force evaluation
  //    std::cout << "(*this)():\n"<<(*this)()<<std::endl;
  //  std::cout << "m_local_representation\n" << m_local_representation << std::endl;
}
GenericOperator Operator::operator*(const Operator& other) const { return GenericOperator(*this, other); }

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
bool Operator::operator==(const Operator& other) const { return not((*this) < other or (other < (*this))); }

bool Operator::operator<(const Operator& other) const {
  constexpr bool debug = false;
  const double tolerance = 1e-12;
  if (debug)
    std::cout << "Operator::operator<: " << this->name() << "\n"
              << (*this)() << ": " << other.name() << "\n"
              << other() << std::endl;
  auto thisu = (*this)();
  auto otheru = other();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (thisu(i, j) < otheru(i, j) - tolerance)
        return true;
      if (thisu(i, j) > otheru(i, j) + tolerance)
        return false;
    }
  }
  return false;
}

std::string GenericOperator::str(const std::string& title, bool coordinate_frame) const {
  std::stringstream result;
  result << "Generic Operator " << title;
  result << Operator::str(title, coordinate_frame);
  result << ", transformation:\n" << this->m_local_representation << std::endl;
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
