import matplotlib as mpl
mpl.use("module://mpl_ascii")

import matplotlib.pyplot as plt

import numpy as np
from tqdm import tqdm

with open("aims.out", "r") as f:
    lines = f.readlines()

energies = []
query = "| Total energy                  :"
for line in tqdm(lines):
    if query in line:
        total_energy_ev = float(line.split()[-2])
        total_energy_ha = float(line.split()[-4])
        energies.append((total_energy_ev, total_energy_ha)))
cycles = len(energies)
energies = np.array(energies)

plt.scatter(range(cycles), energies[:, 0], label="Total DFT energy (eV)")
plt.save
plt.show()

