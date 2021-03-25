#ifndef POINT_CHARGE_SYMMETRY__SYMMETRYOPERATOR_H_
#define POINT_CHARGE_SYMMETRY__SYMMETRYOPERATOR_H_

#include <memory>
#include <array>
#include <string>
#include <ostream>
#include "CoordinateSystem.h"
namespace molpro::point_charge_symmetry {
static CoordinateSystem s_default_coordinate_system;

class Operator {
public:
  using vec = CoordinateSystem::vec;
  using mat = CoordinateSystem::mat;
  virtual ~Operator() = default;
  Operator(const CoordinateSystem& coordinate_system = s_default_coordinate_system)
      : m_coordinate_system(coordinate_system){};
  vec operator()(vec v) const;
  std::array<vec,3> operator_gradient(vec v) const;
  virtual vec operator_local(vec v) const = 0;
  virtual std::string str(const std::string& title) const;
  const std::string& name() const { return m_name; };

protected:
  const CoordinateSystem& m_coordinate_system;
  std::string m_name;
  friend class Group;
};

inline std::ostream& operator<<(std::ostream& os, const Operator& op) {
  os << op.str("");
  return os;
}

class Reflection : public Operator {
protected:
  vec m_normal;

public:
  Reflection(vec normal) : Reflection(s_default_coordinate_system, std::move(normal)) {}
  Reflection(const CoordinateSystem& coordinate_system, vec normal)
      : Operator(coordinate_system), m_normal(normal.normalized()) {
    m_name = "sigma";
    if (m_normal(0) > 0.99)
      m_name += "_yz";
    if (m_normal(1) > 0.99)
      m_name += "_xz";
    if (m_normal(2) > 0.99)
      m_name += "_xy";
  }
  vec operator_local(vec v) const override;
  std::string str(const std::string& title) const override;
  friend class Group;
};

class Rotation : public Operator {
protected:
  vec m_axis;
  int m_order;
  bool m_proper = true;

public:
  Rotation(vec axis, int order = 2, bool proper = true)
      : Rotation(s_default_coordinate_system, std::move(axis), std::move(order), std::move(proper)) {}
  Rotation(const CoordinateSystem& coordinate_system, vec axis, int order = 2, bool proper = true)
      : Operator(coordinate_system), m_axis(axis.normalized()), m_order(std::move(order)), m_proper(std::move(proper)) {
    m_name = (m_proper ? "C" : "S") + std::to_string(m_order);
    if (m_axis(0) > 0.99)
      m_name += "x";
    if (m_axis(1) > 0.99)
      m_name += "y";
    if (m_axis(2) > 0.99)
      m_name += "z";
  }
  vec operator_local(vec v) const override;
  friend class Group;
};

class Inversion : public Operator {
public:
  Inversion(const CoordinateSystem& coordinate_system = s_default_coordinate_system) : Operator(coordinate_system) {
    m_name = "i";
  }
  vec operator_local(vec v) const override;
  friend class Group;
};

class Identity : public Operator {
public:
  Identity(const CoordinateSystem& coordinate_system = s_default_coordinate_system) : Operator(coordinate_system) {
    m_name = "E";
  }
  vec operator_local(vec v) const override;
  friend class Group;
};

} // namespace molpro::point_charge_symmetry

#endif // POINT_CHARGE_SYMMETRY__SYMMETRYOPERATOR_H_
