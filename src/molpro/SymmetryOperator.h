#ifndef POINT_CHARGE_SYMMETRY__SYMMETRYOPERATOR_H_
#define POINT_CHARGE_SYMMETRY__SYMMETRYOPERATOR_H_

#include <molpro/CoordinateSystem.h>
#include <memory>
namespace molpro::point_charge_symmetry {
static CoordinateSystem s_default_coordinate_system;

class SymmetryOperator {
 public:
  using vec = CoordinateSystem::vec;
  using mat = CoordinateSystem::mat;
  virtual ~SymmetryOperator() = default;
  SymmetryOperator(const CoordinateSystem &coordinate_system = s_default_coordinate_system) : m_coordinate_system(coordinate_system) {};
  virtual vec operator()(vec v) const = 0;
  virtual std::string str(const std::string &title) const;
 protected:
  const CoordinateSystem& m_coordinate_system;
  vec global_to_local(vec v) const;
  vec local_to_global(vec v) const;
};

inline std::ostream& operator<<(std::ostream& os, const SymmetryOperator& op) {
  os << op.str("SymmetryOperator");
  return os;
}

class ReflectionPlane : public SymmetryOperator {
 protected:
  CoordinateSystem::vec m_normal;
 public:
  ReflectionPlane(vec normal) : SymmetryOperator(), m_normal(normal.normalized()) {}
  ReflectionPlane(const CoordinateSystem &coordinate_system, vec normal) : SymmetryOperator(
      coordinate_system), m_normal(normal.normalized()) {

  }
  vec operator()(vec v) const override;
  std::string str(const std::string &title) const override;
};

class Axis : public SymmetryOperator {
 protected:
  vec m_axis;
  int m_order;
  bool m_proper = true;
 public:
  Axis(vec axis, int order = 2, bool proper = true)
      : SymmetryOperator(), m_axis(axis.normalized()), m_order(order), m_proper(proper) {}
  Axis(const CoordinateSystem &coordinate_system, vec axis, int order = 2, bool proper = true)
      : SymmetryOperator(coordinate_system), m_axis(axis.normalized()), m_order(order), m_proper(proper) {}
  vec operator()(vec v) const override;
};

class Inversion : public SymmetryOperator {
 public:
  Inversion() : SymmetryOperator() {}
  Inversion(const CoordinateSystem &coordinate_system) : SymmetryOperator(coordinate_system) {}
  vec operator()(vec v) const override;
};
}

#endif //POINT_CHARGE_SYMMETRY__SYMMETRYOPERATOR_H_
