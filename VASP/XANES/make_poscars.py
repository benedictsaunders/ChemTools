from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.build import make_supercell
from ase.build.tools import sort as ase_sort
import numpy as np
from copy import deepcopy as dcp
from tqdm import tqdm
from os import system, environ, path
from data import potentials, XrayNotation
from utils import *
import argparse as ap


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


def read_poscar(inp, outp, P=None):
    at = ase_sort(ase_read(inp))
    if P is not None:
        at = ase_sort(make_supercell(at, P, wrap=False, tol=1e-08))
    ase_write(outp, at, direct=False, wrap=False)

    with open(outp, "r") as f:
        lines = [l.strip() for l in f.readlines()]
        comment = lines[0]
        factor = lines[1]
        lattice = lines[2:5]
        species = lines[5].split()
        counts = [int(l) for l in lines[6].split()]
        coord_type = lines[7]
        coords = lines[8 : sum(counts) + 8]
    poscar = POSCAR(
        lattice=lattice,
        factor=factor,
        species=species,
        counts=counts,
        coord_type=coord_type,
    )
    poscar.sites = []
    ordered_species = []
    for idx, s in enumerate(species):
        ordered_species += [s] * counts[idx]
    for idx, os in enumerate(ordered_species):
        new_site = site(element=os, position=coords[idx], index=idx)
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
    idx = 0
    # Reading the files to POSCAR objects
    prim = read_poscar(inp, "POSCAR_xprimitive")
    supercell = read_poscar(inp, "POSCAR_xsuper", P)

    for site in tqdm(prim.sites, total=len(prim.sites)):
        if site.element == target:
            tqdm.write("Found a target site.")
            idx += 1
            supercellc = dcp(supercell)
            for sidx, supersite in enumerate(supercellc.sites):
                if supersite.position == site.position:
                    print(idx)
                    # supercell_targets[sidx]

                    target_site = supercellc.sites.pop(sidx)
                    supercellc.sites.append(target_site)

                    tidx = supercellc.species.index(target)
                    supercellc.species.append(target)
                    supercellc.counts.append(1)
                    supercellc.counts[tidx] -= 1
                    with cd(f"XANES_{target}_{idx}"):
                        write_poscar(
                            supercellc,
                            name="POSCAR",
                            comment=f" XANES on {target} iteration {idx}",
                        )
                    break
    return (
        len([site.element for site in prim.sites if site.element == target]),
        supercellc.species,
        supercellc.counts,
    )


def make_potcar(order, family="potpaw_PBE", preferred_override=None, useGW=False):
    pppath = environ.get("VASP_PP_PATH") + "/" + family
    pots = []
    for specie in order:
        if preferred_override is not None:
            p = preferred_override[specie]
        else:
            p = potentials.defaults[specie]
        if useGW and path.isfile(f"{pppath}/{p}_GW/POTCAR"):
            p += "_GW"
        pots.append(f"{pppath}/{p}/POTCAR")
    pots_joined = " ".join(pots)
    with open("POTCAR", "w") as potcar:
        _ = runcmd(["cat"] + pots, outlocation=potcar)


def make_incar(
    iters,
    ordering,
    counts,
    target,
    nbands,
    Xtype="K",
    magmoms=None,
    hubbardU=None,
    potpawFamily="potpaw_PBE",
    useGW=False,
    nelect=None,
    P=np.asarray([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
):
    dirs = []
    if not os.path.isfile("INCAR"):
        raise FileNotFoundError("INCAR")
    if nelect is not None:
        detP = np.linalg.det(P)
        nelect = int(nelect) * detP
    for idx in range(1, iters + 1):
        dir = f"XANES_{target}_{idx}"
        dirs.append(dir)
        print(f"Making and populating {dir}.")
        with cd(dir):
            runcmd(["cp", "../INCAR", "."])
            with open("POSCAR", "r") as f:
                lines = [l.strip() for l in f.readlines()]
            species = lines[5].split()
            new_counts = lines[6].split()
            assert len(species) == len(new_counts)
            make_potcar(species, family=potpawFamily, useGW=useGW)
            hubU, hubL, hubJ, mgms = [], [], [], []
            for c, o in zip(counts, ordering):
                if magmoms is not None:
                    try:
                        m = magmoms[o]
                    except:
                        m = 0.0
                    mgms.append(f"{c}*{magmoms[o]}")
                if hubbardU is not None:
                    try:
                        hdict = hubbardU[o]
                    except:
                        hdict = {"l": -1, "U": 0.00, "J": 0.00}
                    hubJ.append(f"{c}*{hdict['J']}")
                    hubU.append(f"{c}*{hdict['U']}")
                    hubL.append(f"{c}*{hdict['l']}")
            with open("INCAR", "a") as f:
                f.write(newline("##### Written by XANESMaker #####"))
                if magmoms is not None:
                    ms = " ".join(mgms)
                    f.write(newline(f"MAGMOM = {ms}"))
                    f.write(newline("ISPIN = 2"))
                if hubbardU is not None:
                    l = " ".join(hubL)
                    U = " ".join(hubU)
                    J = " ".join(hubJ)
                    f.write(newline("LDAU = .TRUE."))
                    f.write(newline("LDAUTYPE = 2"))
                    f.write(newline("LDAUL = " + l))
                    f.write(newline("LDAUU = " + U))
                    f.write(newline("LDAUJ = " + J))
                    f.write(newline("LDAUPRINT = 2"))
                if nelect is not None:
                    f.write(newline(f"NELECT = {nelect}"))
                f.write(newline("ICORELEVEL = 2"))
                f.write(newline(f"CLNT = {len(species)}"))
                f.write(newline(f"CLN = {XrayNotation.edge[Xtype]['n']}"))
                f.write(newline(f"CLL = {XrayNotation.edge[Xtype]['l']}"))
                f.write(newline("CLZ = 1.0"))
                f.write(newline("CH_LSPEC = .TRUE."))
                f.write(newline("CH_SIGMA = 0.5"))
                f.write(newline(f"NBANDS = {nbands}"))
                f.write(newline(""))
    return dirs


def submit(dirs, subcmd, subfile):
    print("Submitting jobs")
    for dir in dirs:
        with cd(dir):
            pwd = os.getcwd()
            _ = runcmd(["cp", f"../{subfile}", "."])
            print(f"Submitting in {pwd}")
            out = runcmd([subcmd, subfile])
            print(out.stdout)


if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument(
        "--supercell", "-r", nargs=3, help="P matrix", default=[2, 2, 2]
    )
    parser.add_argument(
        "--target", "-t", type=str, required=True, help="Target species"
    )
    parser.add_argument("--luj", "-luj", nargs="*", help="Hubbard corrections")
    parser.add_argument("--magmoms", "-mgm", nargs="*", help="Initial magnetic moments")
    parser.add_argument(
        "--input", "-i", help="Input VASP geometry, def.: POSCAR", default="POSCAR"
    )
    parser.add_argument(
        "--simplesiegbahn",
        "-ssieg",
        help="XRay notation of the initial energy level",
        default="K1",
    )
    parser.add_argument(
        "--nelect",
        "-nelect",
        help="Number of electrons in the unit cell",
        default=None,
    )
    parser.add_argument(
        "--bands",
        "-b",
        help="Number of bands to consider in the VASP calculation.",
        default=400,
    )
    parser.add_argument(
        "--sub",
        "-s",
        nargs=2,
        help="Job submission command and file",
        default=["sbatch", "x.sub"],
    )
    parser.add_argument("--gw", "-gw", action="store_true", help="Use GW POTCARs")
    args = parser.parse_args()
    dirs = []
    P = np.array(
        [
            [int(args.supercell[0]), 0, 0],
            [0, int(args.supercell[1]), 0],
            [0, 0, int(args.supercell[2])],
        ]
    )

    input_file = args.input
    target = args.target
    nbands = args.bands
    submission = args.sub

    iters, species_order, species_counts = iterate_supercell_primitive(
        input_file, P=P, target=target
    )
    syms = list(set(species_order))
    magmoms = handle_magmoms(args.magmoms, syms)
    hubbardU = handle_hubbard(args.luj, syms)

    dirs = make_incar(
        iters,
        species_order,
        species_counts,
        target,
        nbands,
        Xtype=args.simplesiegbahn,
        magmoms=magmoms,
        hubbardU=hubbardU,
        useGW=args.gw,
        nelect=args.nelect,
        P=P,
    )

    submit(dirs, submission[0], submission[1])
