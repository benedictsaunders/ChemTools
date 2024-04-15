from ase.io import read, write
atoms = read('periodic/geometry.in')
print(atoms.pbc)
print()
atoms = read('molecule/geometry.in')
print(atoms.pbc)