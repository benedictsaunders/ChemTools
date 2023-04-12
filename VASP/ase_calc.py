"""
ASE_CALC
formerlly VASP 4 USPEX
-- Benedict Saunders, 2022
    A nice little script to use as the vasp runner for USPEX (etc) calculations, which
    allows you to add things like hubbard factors, magetic moments and other such
    properties that are required on a per-atom basis. 
"""

from ase.io.vasp import read_vasp
from ase.calculators.vasp import Vasp
import os

from pprint import pprint

import argparse as ap


def handle_magmoms(atoms, magmoms):
    """
    The magamoms variable should be a list parsed by argparse in the form:
    [...] -mgm Fe 5 Nb 0.6 O 0.6 [...]
    which is then converted to a dictionary:
    d = {
        'Fe': 5.,
        'Nb': 0.6,
        'O': 0.6
        }
    """
    if magmoms == None:
        return atoms
    elements = magmoms[::2]
    values = magmoms[1::2]
    d = dict(zip(elements, values))
    init_mgm = []
    for s in atoms.symbols:
        if s not in elements:
            init_mgm.append(0)
        else:
            init_mgm.append(d[s])
    atoms.set_initial_magnetic_moments(init_mgm)
    return atoms


def handle_hubbard(atoms, luj):
    if luj is None:
        print("Hubbard corrections not set.")
        return None
    labels = ["L", "U", "J"]
    n = 4
    elements = []
    d = {}
    separated = [luj[i: i + n] for i in range(0, len(luj), n)]
    for indiv in separated:
        elements.append(indiv[0])
        d[indiv[0]] = dict(zip(labels, [float(x) for x in indiv[-3:]]))
    for s in atoms.symbols:
        if s not in elements:
            d[s] = dict(
                zip(labels, [-1, 0, 0])
            )  # Negative 1 for no onsite interation added.
    # pprint(d)
    return d


if __name__ == "__main__":

    parser = ap.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        default="POSCAR",
        help="Input filename"
    )
    parser.add_argument(
        "-c",
        "--cores",
        type=int,
        default=128,
        required=False,
        help="Number of cores to run VASP with. Def.: 128",
    )
    parser.add_argument(
        "-vexe",
        "--vasp_executable",
        type=str,
        default="vasp_std",
        required=False,
        help="VASP executable command. Def.: 'vasp_std'",
    )
    parser.add_argument(
        "-pp",
        "--vasp_potentials",
        type=str,
        default="~/APPS/vasp.5.4.1/PPs",
        required=False,
        help="Location of potpaw potentials for VASP",
    )

    parser.add_argument(
        "-mpi",
        "--mpi",
        action="store_true",
        help="Overrides the srun command to mpirun, and changes the syntx of the command accordingly.",
    )

    parser.add_argument(
        "-mgm",
        "--magmoms",
        default=None,
        help="Magnetic moments for a colinear calculation. Eg, 'Fe 5.0 Nb 0.6 O 0.6' Def.: None. If a defined element is not present in the POSCAR, no MAGMOM will be set for it.",
        nargs="*",
        required=False,
    )

    parser.add_argument(
        "-LUJ",
        "--hubbard_correctionLUJ",
        default=None,
        nargs="*",
        required=False,
        help="Hubbard corrections. Usage: <element 1> <L1> <U1> <J1> ... <element n> <Ln> <Un> <Jn>. Def.: None",
    )

    args = parser.parse_args()

    cores = args.cores
    vasp_command = args.vasp_executable
    pp_loc = args.vasp_potentials
    hubbard = args.hubbard_correctionLUJ
    mgms = args.magmoms

    if os.getenv("VASP_PP_PATH") is None:
        os.environ["VASP_PP_PATH"] = f"{pp_loc}"

    if os.getenv("ASE_VASP_COMMAND") is None:
        if not args.mpi:
            os.environ[
                "ASE_VASP_COMMAND"
            ] = f"mpirun -np {cores} {vasp_command}"
        else:
            os.environ["ASE_VASP_COMMAND"] = f"mpirun -np {cores} {vasp_command}"
        print(os.environ["ASE_VASP_COMMAND"])

    atoms = read_vasp("POSCAR")

    # mgms = "Fe 5 Nb 0.6 O 0.6"

    hubbard = handle_hubbard(atoms, args.hubbard_correctionLUJ)

    """
        The below if statement probably can be improved to something that checks for
        ferromagnetism or something similar, rather than being hardcoded
        for Fe. Maybe. Give it a try. I'll still take the credit ;)
    """

    # if any(elem in ["V", "Fe"] for elem in atoms.symbols) and "O" in atoms.symbols:
    #     atoms_tmp = handle_magmoms(atoms, mgms.split(" "))
    #     atoms = atoms_tmp

    """
        If you don't know what's going on here, neither do I. But it works,
        SO DON'T TOUCH IT!
    """
    atoms_tmp = handle_magmoms(atoms, mgms)
    atoms = atoms_tmp

    calc = Vasp()
    calc.read_incar(filename="INCAR")
    calc.set(xc="pbe", ibrion=2, ldau_luj=hubbard)
    atoms.calc = calc
    E = atoms.get_potential_energy()
