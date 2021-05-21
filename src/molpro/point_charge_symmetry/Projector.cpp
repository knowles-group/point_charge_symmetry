#include "Projector.h"
#include <Eigen/core>

namespace molpro::point_charge_symmetry {
Projector::Projector(const Group& group, const Molecule& molecule) : m_n3(molecule.m_atoms.size() * 3) {
  //  std::cout << "group " << group.name() << std::endl;
  auto ginv = double(1) / (group.end() - group.begin());
  Eigen::MatrixXd Q(m_n3, m_n3);
  Q.setZero();
  SymmetryMeasure sm(molecule, group);
  for (const auto& t : group) {
    Eigen::Matrix3d U, identity{Eigen::Matrix3d::Identity()};
    for (int alpha = 0; alpha < 3; alpha++) {
      U.col(alpha) = group.coordinate_system().axes() *
                     t->operator_local((group.coordinate_system().axes().transpose() * identity.col(alpha))).eval();
    }
    //    std::cout << "operator matrix\n" << U << std::endl;
    for (auto A = 0; A < molecule.m_atoms.end() - molecule.m_atoms.begin(); A++) {
      auto C = sm.image_neighbour(A, *t);
      for (int alpha = 0; alpha < 3; alpha++)
        for (int gamma = 0; gamma < 3; gamma++)
          Q(gamma + 3 * C, alpha + 3 * A) += ginv *
                                             //              U(alpha, gamma);
                                             U(gamma, alpha);
    }
  }
  //  std::cout << "Q\n" << m_Q << std::endl;

  Eigen::JacobiSVD<Eigen::MatrixXd> svd;
  svd.compute(Q, Eigen::ComputeFullV);
  svd.setThreshold(0.99);
  //  std::cout << "singular values: " << svd.singularValues().transpose() << std::endl;
  //  std::cout << "rank: " << svd.rank() << std::endl;
  auto vfull = svd.matrixV();
  m_V.assign(vfull.data(), vfull.data() + m_n3 * svd.rank());
//  std::cout << "singular vectors:\n" << Eigen::Map<Eigen::MatrixXd>(m_V.data(), m_n3, svd.rank()).transpose() << std::endl;
}

std::vector<double> Projector::symmetric(std::vector<double> vector) const {
  assert(m_n3 == vector.size());
  auto V = Eigen::Map<const Eigen::MatrixXd>(m_V.data(), m_n3, m_V.size() / m_n3);
  auto result = (V * V.transpose() * Eigen::Map<Eigen::VectorXd>(vector.data(), m_n3)).eval();
  return std::vector<double>(result.eval().data(), result.eval().data() + m_n3);
}

void Projector::remove_symmetric(std::vector<double> vector) const {
  auto symmetric_part = symmetric(vector);
  std::transform(vector.begin(), vector.end(), symmetric_part.begin(), vector.begin(), std::minus<>{});
}
} // namespace molpro::point_charge_symmetry
