#include <array>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iostream>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/point_charge_symmetry/Molecule.h>

using std::cout;
using Rvector = std::array<double, 3>;
TEST(bfgs, quadratic) {
  Rvector x{9e49, 1, -73e8}, g;
  auto solver = molpro::linalg::itsolv::create_Optimize<Rvector, Rvector, Rvector>("BFGS", "max_size_qspace=4");
  int nwork = 1;
  for (int iter = 0; iter < 100; iter++) {
    double value = 0;
    for (int i = 0; i < 3; i++) {
      value += (i + 1) * (x[i] - 1) * (x[i] - 1);
      g[i] = 2 * (i + 1) * (x[i] - 1);
    }
    cout << "x=" << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
    cout << " value " << value << "; gradient=" << g[0] << ", " << g[1] << ", " << g[2] << std::endl;
    nwork = solver->add_value(x, value, g);
    nwork = solver->end_iteration(x, g);
    solver->report();
    if (nwork <= 0)
      break;
  }
}
TEST(bfgs, molecule) {
  molpro::point_charge_symmetry::Molecule water("h2o.xyz");
  cout << water << std::endl;
  auto coc = water.centre_of_charge();
  Rvector x{0, 0, 0}, g;
  auto solver = molpro::linalg::itsolv::create_Optimize<Rvector, Rvector, Rvector>("BFGS", "max_size_qspace=4");
  int nwork = 1;
  for (int iter = 0; iter < 100; iter++) {
    double value = 0;
    std::fill(g.begin(), g.end(), 0);
    for (int pow = 2; pow < 3; pow += 2)
      for (int i = 0; i < 3; i++) {
        value += std::pow(x[i] - coc(i), pow);
        g[i] += pow * std::pow(x[i] - coc(i), pow - 1);
      }
    cout << "x=" << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
    cout << " value " << value << "; gradient=" << g[0] << ", " << g[1] << ", " << g[2] << std::endl;
    nwork = solver->add_value(x, value, g);
    nwork = solver->end_iteration(x, g);
    solver->report();
    if (nwork <= 0)
      break;
  }
}
#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<Rvector, Rvector, Rvector>;