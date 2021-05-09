#!env python3
import sys
import pubchempy
name = sys.argv[1]
filename = name
if (len(sys.argv) > 2):
    filename = sys.argv[2]
results = pubchempy.get_compounds(name,'name',record_type='3d')
for result in results:
    for coord in result.record['coords']:
        for conformer in coord['conformers']:
            print("writing structure to",filename+".xyz")
            with open(filename+".xyz","w") as stream:
                stream.write(str(len(result.elements))+"\n")
                stream.write("PubChem cid="+str(result.cid)+"\n")
                for key, element in enumerate(result.elements):
                    stream.write("{0} {1} {2} {3}\n".format(element,
                    conformer['x'][key],
                    conformer['y'][key],
                    conformer['z'][key])
                    )
            filename += '_'
