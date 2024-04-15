from ase.io import read
from sys import argv
from basis_names import names as basis_dict
import shutil

DEFAULT_SPECIES_LOCATION = "~/software/FHIaims/species_defaults/defaults_2020/"

if __name__ == '__main__':
    if argv[1] in ("-h", "--help"):
        print("Usage: add_basis.py <input_file> <basis_set_level[light, tight, etc.]>")
        exit()
    atoms = read(argv[1])
    species = list(set(atoms.get_atomic_numbers()))
    basis_set_level = argv[2]
    for specie in species:
        basis_file = f"{DEFAULT_SPECIES_LOCATION}{basis_set_level}/{basis_dict[specie]}"
    
