#include "SymmetryMeasureC.h"
#include <molpro/point_charge_symmetry/Group.h>
#include <molpro/point_charge_symmetry/Molecule.h>
#include <molpro/point_charge_symmetry/SymmetryMeasure.h>
#include <string.h>
using molpro::point_charge_symmetry::Group;
using molpro::point_charge_symmetry::Molecule;
double SymmetryMeasureValue(const char *groupname, size_t atoms, const double *coordinates, const double *charges) {
  Eigen::Map<const Eigen::MatrixXd> xyz(coordinates, 3, atoms);
  Eigen::Map<const Eigen::VectorXd> q(charges, atoms);
  auto molecule = Molecule(xyz, q);
  auto group = Group(std::string{groupname});
  molpro::point_charge_symmetry::SymmetryMeasure sm(molecule, group);
  return sm();
}

void SymmetryMeasureOptimiseFrame(const char *groupname, size_t atoms, double *coordinates, const double *charges) {
  Eigen::Map<Eigen::MatrixXd> xyz(coordinates, 3, atoms);
  Eigen::Map<const Eigen::VectorXd> q(charges, atoms);
  const Molecule &molecule0 = Molecule(xyz, q);
  const Group &group = Group(std::string{groupname});
  molpro::point_charge_symmetry::SymmetryMeasure sm(molecule0, group);
  sm.optimise_frame();
  auto molecule = molpro::point_charge_symmetry::molecule_localised(sm.group().coordinate_system(),molecule0);
  for (int i=0; i<atoms; i++)
    xyz.col(i) = molecule.m_atoms[i].position;
}

char *SymmetryMeasureDiscoverGroup(double threshold, size_t atoms, double *coordinates, const double *charges) {
  Eigen::Map<Eigen::MatrixXd> xyz(coordinates, 3, atoms);
  Eigen::Map<const Eigen::VectorXd> q(charges, atoms);
  auto group = molpro::point_charge_symmetry::discover_group(Molecule(xyz,q),threshold);
  return strdup(group.name().c_str());
}
void SymmetryMeasureRefine(const char *groupname, size_t atoms, double *coordinates, const double *charges) {
  Eigen::Map<Eigen::MatrixXd> xyz(coordinates, 3, atoms);
  Eigen::Map<const Eigen::VectorXd> q(charges, atoms);
  molpro::point_charge_symmetry::SymmetryMeasure sm(Molecule(xyz, q), Group(std::string{groupname}));
  auto molecule = sm.refine();
  for (int i=0; i<atoms; i++)
    xyz.col(i) = molecule.m_atoms[i].position;
}
