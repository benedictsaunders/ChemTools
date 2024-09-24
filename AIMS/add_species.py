import os
from sys import argv
from basis_names import names as basis_dict
from basis_names import symbols
import shutil
from copy import deepcopy as dcp

# Default location of species defaults
DEFAULT_SPECIES_LOCATION = f"{os.environ["HOME"]}/software/FHIaims/species_defaults/defaults_2020/"

def concat(infile, outfile, label, append = True):
    """Concatenate the contents of infile to outfile
    Args:
        infile (str): Path to the file to read from
        outfile (str): Path to the file to write to
        append (bool): If True, append to the file, otherwise overwrite it"""
    with open(outfile, 'a' if append else 'w') as out:
        with open(infile, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            newlines = dcp(lines)
            for idx, line in enumerate(lines):
                if line.startswith(f"species"):
                    newlines[idx] = f"species {label}"
                    break
        out.writelines([line + "\n" for line in newlines])

def writeSpeciesToControl(species, level, alltight = False):
    """Write the species defaults to the control.in file
    Args:
        species (list): List of atomic numbers of the species in the system
        level (str): Basis set level to use, must be one of 'light', 'intermediate', 'tight'"""
    
    for s in species:
        label = dcp(s)
        current_level = dcp(level)
        if s.endswith("1"):
            s = s[:-1]
            if not alltight:
                current_level = "light"
        name = basis_dict[symbols[s]]
        print(f"Appending {name} to control.in")
        concat(f"{DEFAULT_SPECIES_LOCATION}{current_level}/{name}", 'control.in', label = label)

def checkExistingSpecies(delete = False):
    """Check if the species defaults are already in the control.in file"""
    existing = []
    with open('control.in', 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.strip().startswith("species"):
                existing.append(line.split()[1])
        print("Existing species:", *[f"\n{spc}" for spc in existing])
    if len(existing) > 0:
        if delete:
            print("Deleting existing species defaults")
            with open('control.in', 'r') as f:
                lines = [line.strip() for line in f.readlines()]

            
        return True
    else:
        return False
    
def read_geometry():
    """Read the geometry from the geometry.in file"""
    species = set()
    with open('geometry.in', 'r') as f:
        for line in f.readlines():
            if line.strip().startswith("atom"):
                species.add(line.strip().split()[4])
    return species
                

if __name__ == '__main__':
    args = argv + [None]
    if argv[1] == None:
        level = "tight"
    else:
        assert argv[1] in ["light", "intermediate", "tight", "alltight"], \
            "Invalid basis set level, must be one of 'light', 'intermediate', 'tight' or 'alltight'"
        level = argv[1]
    species = list(read_geometry())
    if checkExistingSpecies():
        print("Species defaults already in control.in, exiting")
    writeSpeciesToControl(species, level, alltight=(level == "alltight"))


