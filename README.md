## Point Group Discovery and Repair
To use the command-line tool,

``
point_charge_symmetry [--tolerance threshold] [--output filename] [-v] molecule.xyz``

The program reads the molecular geometry in [XYZ](http://openbabel.org/wiki/XYZ_(format)) format, and then tries a succession of candidate point groups. For each one, the symmetry-breaking measure is calculated and minimised with respect to the coordinate frame (origin and orientation) definining the point group. The process stops when the highest order group is found for which the optimum symmetry-breaking measure is less than `threshold` (default 1e-4). 

The geometry is then adjusted with the least-motion displacement that brings the symmetry-breaking measure to zero.

The convenience script `pubchem-xyz.py` does a search by name on [PubChem](https://pubchem.ncbi.nlm.nih.gov), and creates XYZ files for any matching structures.