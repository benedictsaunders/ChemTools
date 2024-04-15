from sys import argv
import numpy as np
from aimstools import KGridFile, GeometryFile, KGridGenerator

def make2DKList(geom, gmax):
    kgrid = KGridFile.KGridFile()
    kgrid.setKGrid(geom=geom, k_grid_type="r_gmp", g_max=gmax, time_reversal=True)
    kgrid.saveToFile("k_list.in")

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