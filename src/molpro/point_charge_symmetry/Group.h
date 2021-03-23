#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_GROUP_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_GROUP_H_
#include <vector>
#include <memory>
#include "Operator.h"

namespace molpro::point_charge_symmetry {
class Group {
 protected:
  CoordinateSystem m_coordinate_system;
  std::vector<std::unique_ptr<Operator>> m_members;
 public:
  Group(const CoordinateSystem &coordinate_system = CoordinateSystem()) : m_coordinate_system(coordinate_system) {}
  void add(Identity op) { m_members.emplace_back(new Identity(m_coordinate_system)); }
  void add(Inversion op) { m_members.emplace_back(new Inversion(m_coordinate_system)); }
  void add(ReflectionPlane op) { m_members.emplace_back(new ReflectionPlane(m_coordinate_system, op.m_normal)); }
  void add(Axis op) { m_members.emplace_back(new Axis(m_coordinate_system, op.m_axis, op.m_order, op.m_proper)); }
  const double *data() const { return m_coordinate_system.data(); }
  double *data() { return m_coordinate_system.data(); }
  using iterator = std::vector<std::unique_ptr<Operator>>::iterator;
  using const_iterator = std::vector<std::unique_ptr<Operator>>::const_iterator;
  iterator begin() { return m_members.begin();}
  const_iterator begin() const { return m_members.begin();}
  iterator end() { return m_members.end();}
  const_iterator end() const { return m_members.end();}
};
}

#endif //POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_GROUP_H_
