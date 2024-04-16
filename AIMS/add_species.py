from ase.io import read
import os
from sys import argv
from basis_names import names as basis_dict
import shutil

# Default location of species defaults
DEFAULT_SPECIES_LOCATION = f"{os.environ["HOME"]}/software/FHIaims/species_defaults/defaults_2020/"

def concat(infile, outfile, append = True):
    """Concatenate the contents of infile to outfile
    Args:
        infile (str): Path to the file to read from
        outfile (str): Path to the file to write to
        append (bool): If True, append to the file, otherwise overwrite it"""
    with open(outfile, 'a' if append else 'w') as out:
        with open(infile, 'r') as f:
            out.write(f.read())

def writeSpeciesToControl(species, level):
    """Write the species defaults to the control.in file
    Args:
        species (list): List of atomic numbers of the species in the system
        level (str): Basis set level to use, must be one of 'light', 'intermediate', 'tight'"""
    for s in species:
        name = basis_dict[s]
        print(f"Appending {name} to control.in")
        concat(f"{DEFAULT_SPECIES_LOCATION}{level}/{name}", 'control.in')

if __name__ == '__main__':
    args = argv + [None]
    if argv[1] == None:
        level = "tight"
    else:
        assert argv[1] in ["light", "intermediate", "tight"], \
            "Invalid basis set level, must be one of 'light', 'intermediate', 'tight'"
        level = argv[1]
    atoms = read('geometry.in')
    species = list(set(atoms.get_atomic_numbers()))
    species.sort()
    writeSpeciesToControl(species, level)


