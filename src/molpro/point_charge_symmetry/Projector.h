#ifndef POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_PROJECTOR_H_
#define POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_PROJECTOR_H_
#include <molpro/point_charge_symmetry/Group.h>
#include <molpro/point_charge_symmetry/Molecule.h>
#include <molpro/point_charge_symmetry/SymmetryMeasure.h>
#include <vector>

namespace molpro::point_charge_symmetry {

class Projector {
public:
  Projector(const Group& group, const Molecule& molecule);
  Projector() = delete;
  std::vector<double> symmetric(std::vector<double> vector) const;
  void remove_symmetric(std::vector<double> vector) const;

protected:
  const size_t m_n3;
  std::vector<double> m_V;
};

} // namespace molpro::point_charge_symmetry
#endif // POINT_CHARGE_SYMMETRY_SRC_MOLPRO_POINT_CHARGE_SYMMETRY_PROJECTOR_H_
