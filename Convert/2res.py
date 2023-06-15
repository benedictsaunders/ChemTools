from ase.io import read as ase_read
import ase.io.res
import argparse as ap
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor as AAA

parser = ap.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-n", "--name", required=True)

args = parser.parse_args()


atoms = ase_read(args.input, index = -1)
struct = AAA.get_structure(atoms)
SGA = SpacegroupAnalyzer(struct)
sg = SGA.get_space_group_symbol()

E = atoms.get_potential_energy()
r = ase.io.res.Res(atoms, energy=E, name=args.name, pressure=0.0, spacegroup = sg)
r.write_file(args.name+".res")
print("Done.")
