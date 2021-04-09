#!env python3
import sys
from chemspipy import ChemSpider
from openbabel import pybel
import os
cs = ChemSpider(os.getenv("CHEMSPIDER_API_KEY"))
name = sys.argv[1]
filename = name
if (len(sys.argv) > 2):
    filename = sys.argv[2]
results = cs.search(name)
for result in results:
    molecule = pybel.readstring("SDF", result.mol_3d)
    while os.path.exists(filename+".xyz"): filename += "_"
    f = open(filename+".xyz","w")
    f.write(molecule.write("XYZ"))
    print("written to "+filename+".xyz")
