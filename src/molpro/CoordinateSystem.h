#ifndef NEAR_SYMMETRY__COORDINATESYSTEM_H_
#define NEAR_SYMMETRY__COORDINATESYSTEM_H_
#include <array>
#include <Eigen/Dense>

namespace molpro::point_charge_symmetry {
/*!
 * @brief Hold the definition of a coordinate system - origin and axes
 */
class CoordinateSystem {
  std::array<double, 6> m_parameters;
 public:
  using vec = Eigen::Vector3d;
  using mat = Eigen::Matrix3d;
  Eigen::Map<vec> origin() { return Eigen::Map<vec>(&m_parameters[0]); }
  Eigen::Map<const vec> origin() const { return Eigen::Map<const vec>(&m_parameters[0]); }
  const mat axes() const;
  Eigen::Map<const vec> axis_generator() const { return Eigen::Map<const vec>(&m_parameters[3]); }
  Eigen::Map<vec> axis_generator() { return Eigen::Map<vec>(&m_parameters[3]); }
  CoordinateSystem(const vec &origin = vec::Zero(), const mat &axes = mat::Identity());
  double *data() { return m_parameters.data(); }
  const double *data() const { return m_parameters.data(); }
};
}

#endif //NEAR_SYMMETRY__COORDINATESYSTEM_H_
