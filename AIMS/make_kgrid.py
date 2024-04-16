from sys import argv
import numpy as np
from aimstools import KGridFile, GeometryFile, KGridGenerator
from copy import deepcopy as dcp

def make2DKList(geometry, gmax):
    kgrid = KGridFile.KGridFile()
    print("Created KGridFile object")

    kgrid.setKGrid(
        geometry,
        k_grid_type="r_gmp",
        g_max=gmax,
        time_reversal=True)
    
    print("Set KGrid")

    kgrid.saveToFile("k_list.in")
    print("Saved KGrid to file")

def print3DGridParams(geom, gmax):
    gmax = gmax * np.pi * 2
    r_lattice = geom.getReciprocalLattice()
    M, N, O = KGridGenerator.getUniformGridParameters3D(r_lattice, gmax)
    print(f"k_grid\t{M}\t{N}\t{O}")

if __name__ == '__main__':
    if argv[1] in ("-h", "--help"):
        print("Usage: make_kgrid.py <geometry.in> <gmax>")
        exit()
    geometry = GeometryFile.GeometryFile()
    geometry.readFromFile(argv[1])
    if argv[2] == "3D":
        print3DGridParams(geometry, float(argv[3]))
    elif argv[2] == "2D":
        make2DKList(geometry, float(argv[3]))