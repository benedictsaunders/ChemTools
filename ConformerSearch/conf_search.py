import os, sys, time, rdkit, argparse
import pickle as pkl
from rdkit import Chem
from rdkit.Chem import rdmolfiles, rdmolops, SanitizeMol, rdDistGeom, AllChem
import matplotlib.pyplot as plt
import pandas as pd

def getXMins(energies, x):
    p = []
    for i in range(x):
        m = min(energies)
        idx = energies.index(m)
        energies[idx] = max(energies) + 1.
        p.append((idx, m))
    print(p)
    return p

def getETKDGconfs(molObj, n=1500, seed=int(time.time()), minimize = True, steps = 1000):
    n = int(n)
    print(f"{n} conformers")
    confenergies = []
    mmffenergies = []
    molObj = Chem.rdmolops.AddHs(molObj)
    molObj = Chem.AddHs(molObj)

    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.randomSeed = int(seed)
    params.numThreads = int(os.environ['OMP_NUM_THREADS'])
    params.useRandomCoords = True
    params.enforceChirality = True

    rmslist1 = []
    rmslist2 = []

    confNum = n
    cids = Chem.rdDistGeom.EmbedMultipleConfs(molObj, confNum, params=params)
    AllChem.AlignMolConformers(molObj, RMSlist=rmslist1)

    # Pickling molecule incase we need to reference it at a later date
    with open('rdkit_mol_raw.pickle', 'wb') as handle:
        pkl.dump(molObj, handle, protocol=pkl.HIGHEST_PROTOCOL)

    for idx, cid in enumerate(cids):
        ff = AllChem.MMFFGetMoleculeForceField(
            molObj, AllChem.MMFFGetMoleculeProperties(molObj, mmffVariant='MMFF94'), confId=cid
        )
        ff.Initialize()
        confenergies.append(ff.CalcEnergy())
    AllChem.AlignMolConformers(molObj, RMSlist=rmslist1)

    res = AllChem.MMFFOptimizeMoleculeConfs(molObj, numThreads = int(os.environ['OMP_NUM_THREADS']))
    AllChem.AlignMolConformers(molObj, RMSlist=rmslist2)

    with open('rdkit_mol_min.pickle', 'wb') as handle:
        pkl.dump(molObj, handle, protocol=pkl.HIGHEST_PROTOCOL)
    for r in res:
        mmffenergies.append(r[1])
    # Get lowest 5 conformers without minimization
    # Get lowest 5 conformer with minimization

    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax1.scatter(cids, confenergies, marker=".")
    ax1.scatter(cids, mmffenergies, marker=".")
    ax1.set_ylabel("Energy $kcal~ mol^{-1}$")
    ax1.legend(["ETKDGv3 output", "MMFF94 minimised"], loc="upper right")

    ax2.scatter(confenergies, mmffenergies, marker=".", color="k")
    ax2.set_ylabel("MMFF94 minimised $kcal~ mol^{-1}$")
    ax2.set_xlabel("ETKDGv3 output $kcal~ mol^{-1}$")
    plt.tight_layout()
    fig.set_size_inches(8, 6)
    plt.savefig("energies.png", transparent=False)
    return molObj, cids, confenergies, mmffenergies

parser = argparse.ArgumentParser(description="Let's do some conformer shenanigans")
parser.add_argument("-s", dest="smiles", action="store", help="SMILES")
parser.add_argument("-n", dest="confs", action="store", help="Must be int.")
args = parser.parse_args()
smiles = args.smiles
n_confs = args.confs
mol = Chem.MolFromSmiles(smiles)
df = pd.DataFrame()
m, cids, confsE, mmffE = getETKDGconfs(mol, n=n_confs)
df['cID'] = cids
df['Raw_conf'] = confsE
df['Min_MMFF'] = mmffE
df.to_csv('conformer_energies.csv')
lowest = getXMins(mmffE, 5)
sorted_min_pairs = sorted(lowest, key=lambda x: x[1])
for idx, pair in enumerate(sorted_min_pairs):
    Chem.rdmolfiles.MolToXYZFile(m, f"{str(idx+1)}_{pair[0]}.xyz", confId=pair[0])
