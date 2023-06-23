from pymatgen.io import cif
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read as ase_read
import numpy as np
import subprocess as sp
import re
from contextlib import contextmanager
import os
import shutil
from glob import glob
import argparse as ap
from templates import FeffInpTemplate

@contextmanager
def cd(path):
    prev = os.getcwd()
    if not os.path.exists(f"{prev}/{path}"):
        os.makedirs(str(path))
    os.chdir(str(path))
    try:
        yield
    finally:
        os.chdir(prev)

def atoms_to_cif(atoms, name="atoms.cif"):
    structure = AseAtomsAdaptor.get_structure(atoms=atoms)
    cw = cif.CifWriter(struct=structure, symprec=0.01)
    cw.write_file(name)

def get_target_indices(name, species):
    targets = []
    with open(name, "r") as f:
        lines = [line.strip() for line in f.readlines()]
    for idx, line in enumerate(lines):
        if "_symmetry_space_group_name_H-M" in line:
            sg = line.split()[-1]
            fmtd_sg = sg.replace("_", "", -1)
            sg_line = idx
        if "_atom_site_occupancy" in line:
            pos_begin = idx + 1
            break
    new_name = name.replace(".cif", "_new.cif")
    lines[sg_line] = "_symmetry_space_group_name_H-M\t " + fmtd_sg
    print(f"Rewriting CIF with formatted spacegroup readable by FEFF as {new_name}")
    with open(new_name, "w") as f:
        f.writelines([line + "\n" for line in lines])

    for target, pos in enumerate(lines[pos_begin:], 1):
        letters = " ".join(re.findall("[a-zA-Z]+", pos)).split()[0]
        if letters == species:
            targets.append(target)
    return targets, new_name

def write_from_template(targets, cif_name, species):
    dirs = []
    with open("feff.inp.template") as template:
        tlines = [tline.strip() for tline in template.readlines()]
    for idx, tline in enumerate(tlines):
        if "TARGET" in tline:
            target_line = idx
        if "CIFLOCATION" in tline:
            cif_line = idx
    for target in targets:
        directory = f"target_{species}_{target}"
        dirs.append(directory)
        with cd(directory):
            tlines[target_line] = f"TARGET  {target}"
            tlines[cif_line] = f"CIF  {cif_name}"
            with open("feff.inp", "w") as inp:
                inp.writelines([tline + "\n" for tline in tlines])
            shutil.copy2(f"../{cif_name}", ".")
    return dirs

def write_from_template_class(title, targets, cif_name):
    dirs = []
    for target in targets:
        directory = f"target_{species}_{target}"
        dirs.append(directory)
        with cd(directory):
            FeffInpTemplate(title, target, cif_name)
            shutil.copy2(f"../{cif_name}", ".")
    return dirs

def run_feff(feff_script, dryrun = False):
    if not dryrun:
        if feff_script[0] != "/":
            print("You must use the absolute path of the FEFF script. Exitting.")
            exit()
        fout = open("feff.out", "w")
        p = sp.Popen(["bash", feff_script], stdout=fout, stderr=fout)
        p.wait()
        fout.close()
        print("<++=== DONE ===++>")
    else:
        print("FEFF DRY RUN")

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--inputs", nargs="+")
    parser.add_argument("-t", "--target", required=True)
    parser.add_argument("--dryrun", action = "store_true")
    args = parser.parse_args()

    target = args.target
    dr = args.dryrun
    files = args.inputs
else:
    files = glob("GLOB HERE")

for f in files:
    atoms = ase_read(f)
    with cd("FEFF_XANES"):
        dname = "".join(f.split("/")[-1].split(".")[:-1])
        with cd(dname):
            species = target
            cifname = f"{dname}.cif"
            atoms_to_cif(atoms, cifname)    
            target_indices, cif_name = get_target_indices(cifname, species)
            #dirs = write_from_template(target_indices, cif_name, species)
            dirs = write_from_template_class(f"{species}_XANES", target_indices, cif_name)
            for d in dirs:
                with cd(d):
                    shutil.copy2(f"../{cif_name}", ".")
                    run_feff(feff_script="/home/chem/msrzvr/scripts/FEFF/feff_mpi.bash", dryrun=dr)