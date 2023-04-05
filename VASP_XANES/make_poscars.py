from sys import argv
from ase.io import read as read_ase
from ase.io import write as write_ase
from pprint import pprint


def make_poscars(sites, input):

    # First, we should normalise the input to the order of things as in ASE
    at = read_ase(input)
    at.sort()
    write_ase(at, "POSCAR_TMP")

    # Reading the lines to designated vars
    with open("POSCAR_TMP", "r") as f:
        lines = [l.strip() for l in f.readlines()]
        comment = lines[0]
        factor = lines[1]
        lattice = lines[2:4]
        species = lines[5].split()
        counts = lines[6].split()
        coord_type = lines[7]
        coords = lines[8:-1]

    coords_of_target = []
    if target not in species:
        print(f"Target species {sites} not found in {input}")
        exit()
    
    d = dict()
    c = 0
    for idx, s in enumerate(species):
        d[s] = (c, counts[idx] + c)
        c =  counts[idx] + c + 1

    pprint(d)

if __name__ == "__main__":
    make_poscars("Fe", argv[1])