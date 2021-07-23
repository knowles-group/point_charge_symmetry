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
    TCLAP::ValueArg<int> repeat("r", "repeat", "Number of executions", false, 1, "integer", cmd);
    TCLAP::ValueArg<double> noise("n", "noise", "Mean noise added to geometries", false, 0, "float", cmd);
    TCLAP::SwitchArg inertial("i", "inertial", "Whether to use inertial axes instead of optimising frame", cmd, false);
    TCLAP::SwitchArg no_frame_optimise("s", "noopt", "Whether to skip frame optimisation", cmd, false);
    TCLAP::SwitchArg project_refine("P", "noproject", "Whether to dskioo projection in structure refinement", cmd,
                                    false);
    TCLAP::ValueArg<double> distance_penalty("d", "distance-penalty", "Penalty on distance in structure refinement",
                                             false, 1e-3, "float", cmd);
    cmd.parse(argc, argv);
    //    std::shared_ptr<molpro::Profiler> prof = molpro::Profiler::single("symmetry_measure");
    //    prof->set_max_depth(1);
    using namespace molpro::point_charge_symmetry;
    molpro::Profiler::single()->set_max_depth(profiler_depth.getValue());
    if (not quiet.getValue() and not given_group.isSet())
      std::cout << "Look for symmetry in " << input_file.getValue() << " with acceptance threshold "
                << tolerance.getValue() << std::endl;
    double total_measure{0};
    for (int count = 0; count < repeat.getValue(); ++count) {
      Molecule molecule(input_file.getValue());
      //      std::cout << "geometry before randomise\n"<<molecule<<std::endl;
      //      std::cout << noise.getValue()<<std::endl;
      molecule.randomise(noise.getValue());
      //      std::cout << "geometry after randomise\n"<<molecule<<std::endl;
      if (not quiet.getValue() and verbose.getValue() > 0)
        std::cout << molecule << std::endl;
      CoordinateSystem cs;
      auto group = given_group.isSet() ? Group(cs, given_group.getValue())
                                       : discover_group(molecule, cs, tolerance.getValue(), verbose.getValue() - 1);
      SymmetryMeasure sm(molecule, group);
      //      std::cout << "symmetry measure before inertial "<<sm()<<std::endl;
      if (not no_frame_optimise.getValue())
        sm.adopt_inertial_axes();
      //      std::cout << "symmetry measure after inertial "<<sm()<<std::endl;
      if (not inertial.getValue() and not no_frame_optimise.getValue())
        sm.refine_frame();
      if (not quiet.getValue() and verbose.getValue() > 3)
        for (size_t i = 0; i < group.size(); i++)
          std::cout << "Symmetry operation " << i << " (" << group[i].name()
                    << ") symmetry-breaking measure = " << sm(i) << std::endl;
      total_measure += sm();
      if (not quiet.getValue() and count == 0) {
        std::cout << (given_group.isSet() ? "Given group " : "Discovered group ") << group.name()
                  << ", symmetry-breaking measure = " << sm() << std::endl;
        if (verbose.getValue() > 1)
          std::cout << group << std::endl;
        if (verbose.getValue() > 1)
          std::cout << cs << std::endl;
      }

      if (not quiet.getValue() and verbose.getValue() > 2)
        std::cout << "Refine geometry " << std::endl;
      auto original_molecule = molecule_localised(cs, molecule);
      //          original_molecule.
      molecule = sm.refine(distance_penalty.getValue(), project_refine.getValue());
      cs = CoordinateSystem();
      if (not quiet.getValue() and count == 0)
        std::cout << "After geometry refinement, symmetry-breaking measure = "
                  << SymmetryMeasure(molecule, Group(group.name()))()
                  << ", cartesian distance from original = " << cartesian_distance(molecule, original_molecule)
                  << ", comparison measure with original = " << distance(molecule, original_molecule).first
                  << std::endl;
      if (not quiet.getValue() and verbose.getValue() > 0)
        std::cout << molecule << std::endl;

      if (output_file.isSet()) {
        if (not quiet.getValue())
          std::cout << "Write geometry to " << output_file.getValue() << std::endl;
        molecule_localised(cs, molecule).write(output_file.getValue());
      }
    }
    if (repeat.getValue() > 1)
      std::cout << "Noise=" << noise.getValue() << ", repeat=" << repeat.getValue()
                << ", mean symmetry measure=" << total_measure / repeat.getValue() << std::endl;
    if (profiler_depth.getValue() > 0)
      std::cout << *molpro::Profiler::single() << std::endl;
  } catch (TCLAP::ArgException& e) // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }
  return 0;
}
