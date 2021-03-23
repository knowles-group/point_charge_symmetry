#include "SymmetryOperator.h"
#include <iostream>
namespace molpro::point_charge_symmetry {
SymmetryOperator::vec SymmetryOperator::global_to_local(vec v) const {
  return (m_coordinate_system.axes().transpose() * (v - m_coordinate_system.origin()));
}
SymmetryOperator::vec SymmetryOperator::local_to_global(vec v) const {
  return (m_coordinate_system.axes() * v + m_coordinate_system.origin());
}

SymmetryOperator::vec ReflectionPlane::operator()(vec v) const {
//  std::cout << "v: "<<v.transpose()<<std::endl;
  auto vlocal = global_to_local(v);
//  std::cout << "vlocal: "<<vlocal.transpose()<<std::endl;
  vlocal -= 2 * m_normal.dot(vlocal) * m_normal;
//  std::cout << "vlocal: "<<vlocal.transpose()<<std::endl;
//  std::cout << "vglobal: "<<local_to_global(vlocal).transpose()<<std::endl;
  return local_to_global(vlocal);
}
SymmetryOperator::vec Axis::operator()(vec v) const {
//  std::cout << "v: "<<v.transpose()<<std::endl;
  auto vlocal = global_to_local(v);
  double angle = (double) 2 * std::acos(double(-1)) / m_order;
  auto aa = Eigen::AngleAxis<double>((double) 2 * std::acos(double(-1)) / m_order, m_axis);
//  std::cout << "vlocal: " << vlocal.transpose() << std::endl;
  vlocal = aa * vlocal;
  if (not m_proper)
    vlocal -= 2 * m_axis.dot(vlocal) * m_axis;
//  std::cout << "vlocal: " << vlocal.transpose() << std::endl;
//  std::cout << "vglobal: " << local_to_global(vlocal).transpose() << std::endl;
  return local_to_global(vlocal);
}

SymmetryOperator::vec Inversion::operator()(vec v) const {
//  std::cout << "v: "<<v.transpose()<<std::endl;
  auto vlocal = global_to_local(v);
//  std::cout << "vlocal: " << vlocal.transpose() << std::endl;
  vlocal = -vlocal;
//  std::cout << "vlocal: " << vlocal.transpose() << std::endl;
//  std::cout << "vglobal: " << local_to_global(vlocal).transpose() << std::endl;
  return local_to_global(vlocal);
}

std::string SymmetryOperator::str(const std::string &title) const {
  std::stringstream result;
  result << "SymmetryOperator";
  result << "\norigin: " << this->m_coordinate_system.origin().transpose();
  result << "\naxes:\n" << this->m_coordinate_system.axes();
  return result.str();
}

std::string ReflectionPlane::str(const std::string &title) const {
  std::stringstream result;
  result << "ReflectionPlane " << title;
  result << SymmetryOperator::str(title);
  result << "\nlocal plane normal: " << this->m_normal.transpose();
  return result.str();
}
}
