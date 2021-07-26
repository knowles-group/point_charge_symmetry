#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_EULER_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_EULER_H_
#include <array>
#include <Eigen/Core>
namespace molpro::point_charge_symmetry {
std::array<double, 3> euler_from_axes(const Eigen::Matrix3d& u);
Eigen::Matrix3d axes_from_euler(const Eigen::Vector3d& euler);
std::array<Eigen::Matrix3d, 3> axes_gradient_from_euler(const Eigen::Vector3d& euler);
}

#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_EULER_H_
