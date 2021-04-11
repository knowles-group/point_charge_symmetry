#include <iostream>
#include <molpro/Profiler.h>
#include <molpro/point_charge_symmetry/CoordinateSystem.h>
#include <molpro/point_charge_symmetry/Molecule.h>
#include <molpro/point_charge_symmetry/SymmetryMeasure.h>
#include <tclap/CmdLine.h>

int main(int argc, char* argv[]) {
  try {
    TCLAP::CmdLine cmd("symmetry_measure", ' ');
    TCLAP::ValueArg<double> tolerance("t", "tolerance", "Maximum symmetry measure for acceptance in group discovery",
                                      false, 1e-4, "tolerance", cmd);
    TCLAP::UnlabeledValueArg<std::string> input_file("file", "Name of xyz file containing structure", true, "", "string",
                                                  cmd);
    TCLAP::MultiSwitchArg verbose("v", "verbose", "show detail", cmd);
    TCLAP::SwitchArg quiet("q", "quiet", "Suppress output", cmd);
    TCLAP::SwitchArg refine("r", "refine", "Refine structure to conform with discovered point group", cmd);
    TCLAP::ValueArg<std::string> output_file("o", "output", "Name of file to contain modified xyz geometry", false, "",
                                        "filename", cmd);
    cmd.parse(argc, argv);
    std::shared_ptr<molpro::Profiler> prof = molpro::Profiler::single("symmetry_measure");
    prof->set_max_depth(1);
    using namespace molpro::point_charge_symmetry;
    if (not quiet.getValue())
      std::cout << "Look for symmetry in " << input_file.getValue() << " with acceptance threshold "
                << tolerance.getValue() << std::endl;
    Molecule molecule(input_file.getValue());
    if (not quiet.getValue() and verbose.getValue() > 0)
      std::cout << molecule << std::endl;
    CoordinateSystem cs;
    auto group = discover_group(molecule, cs, tolerance.getValue(), verbose.getValue()-1);
    if (not quiet.getValue()) {
      std::cout << group.name() << ": " << SymmetryMeasure(molecule, group)() << std::endl;
      if (verbose.getValue() > 1)
        std::cout << group << std::endl;
    }

    if (refine.isSet()) {
      if (not quiet.getValue())
        std::cout << "Refine geometry " << std::endl;
      molecule = SymmetryMeasure(molecule, group).refine(3);
      cs = CoordinateSystem();
      if (not quiet.getValue())
        std::cout << "Refined symmetry measure = " << SymmetryMeasure(molecule, Group(group.name()))() << std::endl;
      std::cout <<molecule<<std::endl;
    }

    if (output_file.isSet()) {
      if (not quiet.getValue())
        std::cout << "Write geometry to " << output_file.getValue() << std::endl;
      molecule_localised(cs,molecule).write(output_file.getValue());
    }

  } catch (TCLAP::ArgException& e) // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }
  return 0;
}
