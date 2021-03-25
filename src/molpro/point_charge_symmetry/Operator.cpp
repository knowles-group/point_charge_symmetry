#include "Operator.h"
#include <memory>
#include <sstream>
#include <vector>
//#include <iostream>
namespace molpro::point_charge_symmetry {

Operator::vec Operator::operator()(vec v) const {
  return m_coordinate_system.origin() +
         m_coordinate_system.axes() *
             operator_local((m_coordinate_system.axes().transpose() * (v - m_coordinate_system.origin()))).eval();
}
std::array<Operator::vec, 6> Operator::operator_gradient(vec v, int numerical, double step) const {
  std::array<Operator::vec, 6> result;
  if (numerical > 0) {
    CoordinateSystem coordinate_system(m_coordinate_system);
    for (int i = 0; i < 6; i++) {
      std::vector<std::unique_ptr<Operator>> points;
      std::vector<CoordinateSystem> coordinate_systems;
      std::vector<vec> transformed_points;
      for (int displacement = -numerical; displacement <= numerical; displacement++) {
        coordinate_systems.emplace_back(m_coordinate_system);
        coordinate_systems.back().data()[i] += step * displacement;
        points.emplace_back(this->clone(coordinate_systems.back()));
        transformed_points.push_back((*points.back())(v).eval());
      }
      if (numerical == 1)
//        result[i] = ((*points[2])(v) - (*points[0])(v)).eval() / (2 * step);
      result[i] = (transformed_points[2] - transformed_points[0]) / (2 * step);
      else if (numerical == 2)
        result[i] = (transformed_points[0]-8*transformed_points[1]+8*transformed_points[3]-transformed_points[4]) / (12 * step);
      else
        throw std::logic_error("Incorrect differentiation order");
    }
//    std::cout << "result ";
//    for (int i=0; i<6; i++) std::cout <<" "<<result[i];std::cout<<std::endl;
    return result;
  }
  return result;
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
