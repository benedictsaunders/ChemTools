from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.build import make_supercell
from ase.build.tools import sort as ase_sort
from ase import Atoms
import numpy as np
from sys import argv
from pprint import pprint
from copy import deepcopy as dcp
from tqdm import tqdm

def newline(s):
    if type(s) == list:
        l = [str(x) for x in s]
        s = " ".join(l)
    return str(s) + "\n"

class site:
    def __init__(self, element, position, index) -> None:
        self.element = element
        self.position = position
        self.index = index

class POSCAR:
    def __init__(self, lattice, factor, species, counts, coord_type) -> None:
        self.lattice = lattice
        self.factor = factor
        self.sites = []
        self.species = species
        self.counts = counts
        self.coord_type = coord_type

def read_poscar(inp, outp, P = None):
    at = ase_sort(ase_read(inp))
    if P is not None:
        at = ase_sort(make_supercell(at, P))
    ase_write(outp, at)

    with open(outp, "r") as f:
            lines = [l.strip() for l in f.readlines()]
            comment = lines[0]
            factor = lines[1]
            lattice = lines[2:5]
            species = lines[5].split()
            counts = [int(l) for l in lines[6].split()]
            coord_type = lines[7]
            coords = lines[8:sum(counts)+8]
    poscar = POSCAR(
        lattice=lattice,
        factor=factor,
        species=species,
        counts = counts,
        coord_type=coord_type
    )
    poscar.sites = []
    ordered_species = []
    for idx, s in enumerate(species):
        ordered_species += [s] * counts[idx]
    for idx, os in enumerate(ordered_species):
        new_site = site(
            element=os,
            position=coords[idx],
            index=idx
        )
        poscar.sites.append(new_site)
        del new_site

    return poscar

def write_poscar(poscar, name="POSCAR", comment=""):
    lines = []
    lines.append(newline(comment))
    lines.append(newline(poscar.factor))
    lines.append(newline(poscar.lattice[0]))
    lines.append(newline(poscar.lattice[1]))
    lines.append(newline(poscar.lattice[2]))
    lines.append(newline(poscar.species))
    lines.append(newline(poscar.counts))
    lines.append(newline(poscar.coord_type))
    for site in poscar.sites:
        lines.append(newline(site.position))

    with open(name, "w") as f:
        f.writelines(lines)


def iterate_supercell_primitive(inp, P, target):
     
    prim = read_poscar(inp, "POSCAR_xprimitive")
    supercell = read_poscar(inp, "POSCAR_xsuper", P)

    supercell_targets = []

    for idx, site in tqdm(enumerate(prim.sites, 1), total=len(prim.sites)):
        if site.element == target:
            supercellc = dcp(supercell)
            for sidx, supersite in enumerate(supercellc.sites):
                if supersite.position == site.position:
                    #supercell_targets[sidx] 

                    target_site = supercellc.sites.pop(sidx)
                    supercellc.sites.append(target_site)

                    tidx = supercellc.species.index(target)
                    supercellc.species.append(target)
                    supercellc.counts.append(1)
                    supercellc.counts[tidx] -= 1
                    write_poscar(supercellc, name=f"POSCAR_i{idx}", comment=f" XANES on {target} iteration {idx}")


        


if __name__ == "__main__":
    P = np.array([
        [2,0,0],[0,2,0],[0,0,2]
    ])
    c = iterate_supercell_primitive(argv[1], P=P, target=argv[2])

