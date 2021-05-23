#include "SymmetryMeasureC.h"
#include <molpro/point_charge_symmetry/Group.h>
#include <molpro/point_charge_symmetry/Molecule.h>
#include <molpro/point_charge_symmetry/SymmetryMeasure.h>
#include <string.h>
using molpro::point_charge_symmetry::Group;
using molpro::point_charge_symmetry::Molecule;
using molpro::point_charge_symmetry::CoordinateSystem;
double SymmetryMeasureValue(const char *groupname, size_t atoms, const double *coordinates, const double *charges) {
  Eigen::Map<const Eigen::MatrixXd> xyz(coordinates, 3, atoms);
  Eigen::Map<const Eigen::VectorXd> q(charges, atoms);
  auto molecule = Molecule(xyz, q);
//  std::cout << "in SymmetryMeasureValue, molecule="<<molecule<<std::endl;
  CoordinateSystem cs;
  auto group = Group(cs, std::string{groupname});
  molpro::point_charge_symmetry::SymmetryMeasure sm(molecule, group);
  return sm();
}

void SymmetryMeasureOptimiseFrame(const char *groupname, size_t atoms, double *coordinates, const double *charges) {
  Eigen::Map<Eigen::MatrixXd> xyz(coordinates, 3, atoms);
  Eigen::Map<const Eigen::VectorXd> q(charges, atoms);
  const Molecule &molecule0 = Molecule(xyz, q);
  CoordinateSystem cs;
  const auto group = Group(cs, std::string{groupname});
  molpro::point_charge_symmetry::SymmetryMeasure sm(molecule0, group);
  sm.adopt_inertial_axes();
  sm.refine_frame();
//  std::cout << "in SymmetryMeasureOptimiseFrame after refine_frame() "<<sm()<<std::endl;
  auto molecule = molpro::point_charge_symmetry::molecule_localised(sm.group().coordinate_system(),molecule0);
//  CoordinateSystem cs2;
//  const auto group2 = Group(cs2, std::string{groupname});
//  molpro::point_charge_symmetry::SymmetryMeasure sm2(molecule, group2);
//  std::cout << "in SymmetryMeasureOptimiseFrame after molecule_localised() "<<sm2()<<std::endl;
//  std::cout << "in SymmetryMeasureOptimiseFrame, molecule="<<molecule<<std::endl;
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
  CoordinateSystem cs;
  const auto group = Group(cs, std::string{groupname});
  auto molecule1 = Molecule(xyz, q);
  molpro::point_charge_symmetry::SymmetryMeasure sm(molecule1, group);
  auto molecule = sm.refine();
  for (int i=0; i<atoms; i++)
    xyz.col(i) = molecule.m_atoms[i].position;
//  CoordinateSystem cs2;
//  const auto group2 = Group(cs2, std::string{groupname});
//  molpro::point_charge_symmetry::SymmetryMeasure sm2(molecule, group2);
//  std::cout << "afer refinement "<<sm2()<<", "<<SymmetryMeasureValue(groupname,atoms,coordinates,charges)<<std::endl;
}
