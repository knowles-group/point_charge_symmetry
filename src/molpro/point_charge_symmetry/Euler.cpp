#include "Euler.h"
//#include <iostream>
namespace molpro::point_charge_symmetry {
std::array<double, 3> euler_from_axes(const Eigen::Matrix3d& u) {
  //  std::cout << "euler_from_axes\n" << u << std::endl;
  const auto pi = std::acos(double(-1));
  double alpha, beta, gamma;
  beta = std::acos(std::min(double(1), std::max(double(-1), u(2, 2))));
  if (beta < 1e-6)
    beta = std::asin(std::sqrt(std::pow(u(0, 2), 2) + std::pow(u(1, 2), 2)));
  if (beta < 1e-12) {
    alpha = std::atan2(u(1, 0), u(0, 0));
    gamma = 0;
  } else if (u(0, 0) > 1 - 1e-10) {
    alpha = pi / 2;
    gamma = -pi / 2;
  } else {
    alpha = std::atan2(u(1, 2), u(0, 2));
    gamma = std::atan2(u(2, 1), -u(2, 0));
    if (gamma >= pi)
      gamma -= pi;
  }
  //  std::cout << "Euler angles: " << alpha << ", " << beta << ", " << gamma << std::endl;
  return {alpha, beta, gamma};
}

std::array<Eigen::Matrix3d, 3> axes_gradient_from_euler(const Eigen::Vector3d& euler) {
  std::array<Eigen::Matrix3d, 3> result;
  auto c1 = std::cos(euler[0]);
  auto s1 = std::sin(euler[0]);
  auto c2 = std::cos(euler[1]);
  auto s2 = std::sin(euler[1]);
  auto c3 = std::cos(euler[2]);
  auto s3 = std::sin(euler[2]);
  result[0] << -c2 * c3 * s1 - c1 * s3, -c1 * c3 + c2 * s1 * s3, -s1 * s2, c1 * c2 * c3 - s1 * s3,
      -c3 * s1 - c1 * c2 * s3, c1 * s2, 0, 0, 0;
  result[1] << -c1 * c3 * s2, c1 * s2 * s3, c1 * c2, -c3 * s1 * s2, s1 * s2 * s3, c2 * s1, -c2 * c3, c2 * s3, -s2;
  result[2] << -c3 * s1 - c1 * c2 * s3, -c1 * c2 * c3 + s1 * s3, 0, c1 * c3 - c2 * s1 * s3, -c2 * c3 * s1 - c1 * s3, 0,
      s2 * s3, c3 * s2, 0;
  return result;
}

Eigen::Matrix3d axes_from_euler(const Eigen::Vector3d& euler) {
  auto c1 = std::cos(euler[0]);
  auto s1 = std::sin(euler[0]);
  auto c2 = std::cos(euler[1]);
  auto s2 = std::sin(euler[1]);
  auto c3 = std::cos(euler[2]);
  auto s3 = std::sin(euler[2]);
  Eigen::Matrix3d result;
  result << c1 * c2 * c3 - s1 * s3, -c3 * s1 - c1 * c2 * s3, c1 * s2, c2 * c3 * s1 + c1 * s3, c1 * c3 - c2 * s1 * s3,
      s1 * s2, -c3 * s2, s2 * s3, c2;
  if (false) {
    if (result(0, 0) < 0 and result(1, 1) < 0) {
      result.col(0) = -result.col(0);
      result.col(1) = -result.col(1);
    }
    if (result(0, 0) < 0 and result(2, 2) < 0) {
      result.col(0) = -result.col(0);
      result.col(2) = -result.col(2);
    }
    if (result(1, 1) < 0 and result(2, 2) < 0) {
      result.col(1) = -result.col(1);
      result.col(2) = -result.col(2);
    }
  }
  return result;
}

} // namespace molpro::point_charge_symmetry
