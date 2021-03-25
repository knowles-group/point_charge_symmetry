#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_GROUP_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_GROUP_H_
#include "Operator.h"
#include <memory>
#include <string>
#include <vector>

namespace molpro::point_charge_symmetry {
class Group {
protected:
  CoordinateSystem m_coordinate_system;
  std::vector<std::unique_ptr<Operator>> m_members;
  const std::string m_name;

public:
  Group(const CoordinateSystem &coordinate_system = CoordinateSystem(), std::string name = "")
      : m_coordinate_system(coordinate_system), m_name(std::move(name)) {}
  Group(std::string name) : m_coordinate_system(CoordinateSystem()), m_name(std::move(name)) {}
  std::string name() const { return m_name; }
  void add(Identity op) { m_members.emplace_back(new Identity(m_coordinate_system)); }
  void add(Inversion op) { m_members.emplace_back(new Inversion(m_coordinate_system)); }
  void add(Reflection op) { m_members.emplace_back(new Reflection(m_coordinate_system, op.m_normal)); }
  void add(Rotation op) { m_members.emplace_back(new Rotation(m_coordinate_system, op.m_axis, op.m_order, op.m_proper)); }
//  const double *data() const { return m_coordinate_system.data(); }
//  double *data() { return m_coordinate_system.data(); }
  const CoordinateSystem coordinate_system() const { return m_coordinate_system; }
  using iterator = std::vector<std::unique_ptr<Operator>>::iterator;
  using const_iterator = std::vector<std::unique_ptr<Operator>>::const_iterator;
  iterator begin() { return m_members.begin(); }
  const_iterator begin() const { return m_members.begin(); }
  iterator end() { return m_members.end(); }
  const_iterator end() const { return m_members.end(); }
  Operator &operator[](size_t index) { return *(m_members[index]); }
  const Operator &operator[](size_t index) const { return *(m_members[index]); }
};
} // namespace molpro::point_charge_symmetry

#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_GROUP_H_
