#ifndef NEAR_SYMMETRY__COORDINATESYSTEM_H_
#define NEAR_SYMMETRY__COORDINATESYSTEM_H_
#include <Eigen/Dense>
#include <array>
#include <ostream>

namespace molpro::point_charge_symmetry {
/*!
 * @brief Hold the definition of a coordinate system - origin and axes
 */
class CoordinateSystem {
public:
  using parameters_t = std::array<double, 6>;
  // protected:
  mutable parameters_t m_parameters;

public:
  using vec = Eigen::Vector3d;
  using mat = Eigen::Matrix3d;
  Eigen::Map<vec> origin() { return Eigen::Map<vec>(&m_parameters[0]); }
  Eigen::Map<const vec> origin() const { return Eigen::Map<const vec>(&m_parameters[0]); }
  const mat axes() const;
  const std::array<mat, 3> axes_gradient(int displacements = 2, double step = 1e-4) const;
  Eigen::Map<const vec> axis_generator() const { return Eigen::Map<const vec>(&m_parameters[3]); }
  Eigen::Map<vec> axis_generator() { return Eigen::Map<vec>(&m_parameters[3]); }
  CoordinateSystem(const vec& origin = vec::Zero(), const mat& axes = mat::Identity());
  double* data() { return m_parameters.data(); }
  const double* data() const { return m_parameters.data(); }
  std::string str() const;
  vec to_local(const vec& source);
  vec to_global(const vec& source);
  void cycle_axes();
  void rot90(int axis);
};

inline std::ostream& operator<<(std::ostream& os, const CoordinateSystem& op) {
  os << op.str();
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const CoordinateSystem::parameters_t& p) {
  for (const auto& pp : p)
    os << " " << pp;
  return os;
}

} // namespace molpro::point_charge_symmetry

#endif // NEAR_SYMMETRY__COORDINATESYSTEM_H_
