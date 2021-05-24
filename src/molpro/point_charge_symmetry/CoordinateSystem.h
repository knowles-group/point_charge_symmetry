#ifndef NEAR_SYMMETRY__COORDINATESYSTEM_H_
#define NEAR_SYMMETRY__COORDINATESYSTEM_H_
#include <Eigen/Dense>
#include <array>
#include <ostream>

namespace molpro::point_charge_symmetry {
enum class RotationParameterType { Log, Euler, Quaternion, TanEuler };
/*!
 * @brief Hold the definition of a coordinate system - origin and axes
 */
class CoordinateSystem {
public:
  using parameters_t = std::array<double, 6>;
  // protected:
  mutable parameters_t m_parameters;
  const RotationParameterType m_rotation_parameter_type = RotationParameterType::Euler;
  const double m_rotation_parameter_scale = double(0*1e-11)+1;

protected:
  mutable bool m_axis_permutation_rot90_next = false;

public:
  using vec = Eigen::Vector3d;
  using mat = Eigen::Matrix3d;
  CoordinateSystem& operator=(const CoordinateSystem& source) ;
  Eigen::Map<vec> origin() { return Eigen::Map<vec>(&m_parameters[0]); }
  Eigen::Map<const vec> origin() const { return Eigen::Map<const vec>(&m_parameters[0]); }
  const mat axes() const;
  const std::array<mat, 3> axes_gradient(int displacements = 0, double step = 1e-4) const;
  Eigen::Map<const vec> axis_generator() const { return Eigen::Map<const vec>(&m_parameters[3]); }
  Eigen::Map<vec> axis_generator() { return Eigen::Map<vec>(&m_parameters[3]); }
  CoordinateSystem(const RotationParameterType rotation_parameter_type = RotationParameterType::Euler,
                   const vec& origin = vec::Zero(), const mat& axes = mat::Identity());
  void from_axes(const mat& axes = mat::Identity()) const;
  double* data() { return m_parameters.data(); }
  const double* data() const { return m_parameters.data(); }
  std::string str() const;
  vec to_local(const vec& source) const;
  vec to_global(const vec& source) const;
  void cycle_axes() const;
  std::array<std::array<double,2>,3> rotation_generator_ranges() const;
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
