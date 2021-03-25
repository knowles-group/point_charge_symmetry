#include "Operator.h"
#include <sstream>
namespace molpro::point_charge_symmetry {

Operator::vec Operator::operator()(vec v) const {
  return m_coordinate_system.origin() +
         m_coordinate_system.axes() *
             operator_local((m_coordinate_system.axes().transpose() * (v - m_coordinate_system.origin())));
}
Operator::vec Reflection::operator_local(vec v) const {
  v -= 2 * m_normal.dot(v) * m_normal;
  return v;
}
Operator::vec Rotation::operator_local(vec v) const {
  double angle = (double)2 * std::acos(double(-1)) / m_order;
  auto aa = Eigen::AngleAxis<double>((double)2 * std::acos(double(-1)) / m_order, m_axis);
  v = aa * v;
  if (not m_proper)
    v -= 2 * m_axis.dot(v) * m_axis;
  return v;
}

Operator::vec Inversion::operator_local(vec v) const { return -v; }

Operator::vec Identity::operator_local(vec v) const { return v; }

std::string Operator::str(const std::string& title) const {
  std::stringstream result;
  result << "SymmetryOperator " + m_name;
  if (!title.empty())
    result << " " << title;
  result << "\norigin: " << this->m_coordinate_system.origin().transpose();
  result << "\naxes:\n" << this->m_coordinate_system.axes();
  return result.str();
}

std::string Reflection::str(const std::string& title) const {
  std::stringstream result;
  result << "ReflectionPlane " << title;
  result << Operator::str(title);
  result << "\nlocal plane normal: " << this->m_normal.transpose();
  return result.str();
}
} // namespace molpro::point_charge_symmetry
