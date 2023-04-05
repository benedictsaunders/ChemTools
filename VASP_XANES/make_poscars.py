from sys import argv
from ase.io import read as read_ase
from ase.io import write as write_ase
from ase.build.tools import sort as ase_sort
from pprint import pprint
from copy import deepcopy as dcp

def newline(s):
    return s+ "\n"

def write_poscar(name, comment, factor, lattice, species, counts, coord_type, positions):
    with open(name, "w") as f:
        f.writelines([str(s) for s in [
            newline(comment),
            newline(factor),
            newline(lattice[0]),
            newline(lattice[1]),
            newline(lattice[2]),
            newline(species),
            newline(counts),
            newline(coord_type),
        ]] + [newline(p) for p in positions])
    return 0

def make_poscars(target, inp):

    # First, we should normalise the inp to the order of things as in ASE
    at = read_ase(inp)
    at = ase_sort(at)
    write_ase("POSCAR_TMP", at)

    # Reading the lines to designated vars
    with open("POSCAR_TMP", "r") as f:
        lines = [l.strip() for l in f.readlines()]
        comment = lines[0]
        factor = lines[1]
        lattice = lines[2:5]
        species = lines[5].split()
        counts = [int(l) for l in lines[6].split()]
        coord_type = lines[7]
        coords = lines[8:-1]

    # Getting the indexed positions of each specie in the POSCAR
    # and reordering the counts and species.
    coords_of_target = []
    if target not in species:
        print(f"Target species {target} not found in {inp}")
        exit()
    
    d = dict()
    c = 0
    for idx, s in enumerate(species):
        if s == target:
            new_species = species + [target]
            target_pos = idx
        try:
            d[s] = d[s] + list(range(c, counts[idx] + c))
        except:
            d[s] = list(range(c, counts[idx] + c))
        c =  counts[idx] + c + 1

    counts.append(1)
    counts[target_pos] -= 1

    # Iterating the target positions to move repsective coordinates to end of list
    print(d[target])
    for idx, atom_idx in enumerate(d[target]):
        new_coords = dcp(coords)
        new_coords.append(new_coords.pop(atom_idx))
        _ = write_poscar(
            f"{inp}_{idx+1}",
            comment + f" XANES on {target} iteration {idx+1}",
            factor,
            lattice,
            " ".join(new_species),
            " ".join(str(c) for c in counts),
            coord_type,
            new_coords,
        )

if __name__ == "__main__":
    make_poscars("Fe", argv[1])