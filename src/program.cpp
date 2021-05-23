#include <iostream>
//#include <molpro/Profiler.h>
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
    TCLAP::UnlabeledValueArg<std::string> input_file("file", "Name of xyz file containing structure", true, "",
                                                     "string", cmd);
    TCLAP::MultiSwitchArg verbose("v", "verbose", "show detail", cmd);
    TCLAP::SwitchArg quiet("q", "quiet", "Suppress output", cmd);
    TCLAP::ValueArg<std::string> given_group("g", "group",
                                             "Name of point group to test; if omitted, the group is instead discovered",
                                             false, "", "Mixed case Schoenflies notation", cmd);
    TCLAP::ValueArg<std::string> output_file("o", "output", "Name of file to contain modified xyz geometry", false, "",
                                             "filename", cmd);
    TCLAP::ValueArg<int> profiler_depth("p", "profiler", "Maximum depth of execution time profiling", false, -1,
                                        "integer", cmd);
    cmd.parse(argc, argv);
    //    std::shared_ptr<molpro::Profiler> prof = molpro::Profiler::single("symmetry_measure");
    //    prof->set_max_depth(1);
    using namespace molpro::point_charge_symmetry;
    molpro::Profiler::single()->set_max_depth(profiler_depth.getValue());
    if (not quiet.getValue() and not given_group.isSet())
      std::cout << "Look for symmetry in " << input_file.getValue() << " with acceptance threshold "
                << tolerance.getValue() << std::endl;
    Molecule molecule(input_file.getValue());
    if (not quiet.getValue() and verbose.getValue() > 0)
      std::cout << molecule << std::endl;
    CoordinateSystem cs;
    auto group = given_group.isSet() ? Group(cs, given_group.getValue())
                                     : discover_group(molecule, cs, tolerance.getValue(), verbose.getValue() - 1);
    SymmetryMeasure sm(molecule, group);
    sm.refine_frame();
    if (not quiet.getValue() and verbose.getValue() > 3)
      for (int i = 0; i < group.end() - group.begin(); i++)
        std::cout << "Symmetry operation " << i << " (" << group[i].name() << ") symmetry-breaking measure = " << sm(i)
                  << std::endl;
    if (not quiet.getValue()) {
      std::cout << (given_group.isSet() ? "Given group " : "Discovered group ") << group.name()
                << ", symmetry-breaking measure = " << sm() << std::endl;
      if (verbose.getValue() > 1)
        std::cout << group << std::endl;
      if (verbose.getValue() > 1)
        std::cout << cs << std::endl;
    }

    if (not quiet.getValue() and verbose.getValue() > 2)
      std::cout << "Refine geometry " << std::endl;
    molecule = sm.refine(3);
    cs = CoordinateSystem();
    if (not quiet.getValue())
      std::cout << "After geometry refinement, symmetry-breaking measure = "
                << SymmetryMeasure(molecule, Group(group.name()))() << std::endl;
    if (not quiet.getValue() and verbose.getValue() > 0)
      std::cout << molecule << std::endl;

    if (output_file.isSet()) {
      if (not quiet.getValue())
        std::cout << "Write geometry to " << output_file.getValue() << std::endl;
      molecule_localised(cs, molecule).write(output_file.getValue());
    }
    if (profiler_depth.getValue() > 0)
      std::cout << *molpro::Profiler::single() << std::endl;
  } catch (TCLAP::ArgException& e) // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }
  return 0;
}
