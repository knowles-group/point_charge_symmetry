#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_GROUP_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_GROUP_H_
#include "Operator.h"
#include <memory>
#include <set>
#include <string>

namespace molpro::point_charge_symmetry {
template <class T>
class lesstarget {
public:
  bool operator()(const T& A, const T& B) const { return (*A) < (*B); }
};
static CoordinateSystem s_group_default_coordinate_system;
class Group {
protected:
  CoordinateSystem& m_coordinate_system;
  std::set<std::unique_ptr<Operator>
      , lesstarget<std::unique_ptr<Operator>>
           > m_members;
  std::string m_name;

public:
  Group();
  Group(CoordinateSystem& coordinate_system);
  Group(CoordinateSystem& coordinate_system, std::string name, bool generators_only = false);
  Group(const std::string& name, bool generators_only = false);
  //  Group(std::string name);
  Group(CoordinateSystem& coordinate_system, const Group& source);
  std::string name() const { return m_name; }
  std::string& name() { return m_name; }
  void add(GenericOperator op) { m_members.emplace(new GenericOperator(op)); }
  void add(Identity op) { m_members.emplace(new Identity(m_coordinate_system)); }
  void add(Inversion op) { m_members.emplace(new Inversion(m_coordinate_system)); }
  void add(Reflection op) { m_members.emplace(new Reflection(m_coordinate_system, op.m_normal)); }
  void add(const Rotation& op) {
    m_members.emplace(new Rotation(m_coordinate_system, op.m_axis, op.m_order, op.m_proper, op.m_count));
  }
  void clear() { m_members.clear(); }
  //  const double *data() const { return m_coordinate_system.data(); }
  //  double *data() { return m_coordinate_system.data(); }
  const CoordinateSystem& coordinate_system() const { return m_coordinate_system; }
  CoordinateSystem::parameters_t& coordinate_system_parameters() const { return m_coordinate_system.m_parameters; }
  using iterator = std::set<std::unique_ptr<Operator>>::iterator;
  using const_iterator = std::set<std::unique_ptr<Operator>>::const_iterator;
  iterator begin() { return m_members.begin(); }
  const_iterator begin() const { return m_members.begin(); }
  iterator end() { return m_members.end(); }
  const_iterator end() const { return m_members.end(); }
  size_t size() const { return m_members.size(); }
  Rotation highest_rotation(bool proper = true, size_t index = 0) const;
  //  Operator& operator[](size_t index) { return *(m_members[index]); }
  const Operator& operator[](size_t index) const;
  const_iterator find(const Operator& key) const {
    for (auto member = m_members.begin(); member != m_members.end(); ++member) {
      if (*(*member) == key)
        return member;
    }
    throw std::runtime_error("bad key");
  }
};
inline std::ostream& operator<<(std::ostream& os, const Group& g) {
  os << "Group " << g.name() << ", order=" << g.size() << "\n" << g.coordinate_system();
  for (const auto& el : g)
    os << "\nMember: " << *el;
  return os;
}

Group generate(const Group& generator);

} // namespace molpro::point_charge_symmetry

#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_GROUP_H_
