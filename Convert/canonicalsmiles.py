import sqlite3, string, argparse
import numpy as np
import pandas as pd
from rdkit import Chem

def toCanonical(s):
    molObj = Chem.MolFromSmiles(s)
    Chem.SanitizeMol(molObj)
    return Chem.MolToSmiles(molObj, canonical=True)

parser = argparse.ArgumentParser()
parser.add_argument('-d', dest='database', default='data.db')
parser.add_argument('-t', dest='table')
targs = parser.parse_args()

connection = sqlite3.connect(targs.database)
table = targs.table

df = pd.read_sql_query(f"select * from {table}", connection)
cols = list(df.columns)
print(cols)

cancols = []
mons = []

for col in cols:
    if col in string.ascii_uppercase:
        mons.append(col)
        canonical = []
        smiles = df[col]
        for s in smiles:
            canonical.append(toCanonical(s))
        cancols.append(canonical)

loc = 2
for idx, c in enumerate(cancols, loc):
    df.insert(idx, f"Canonical_{mons[idx-loc]}", c)
df.to_sql(f"{targs.table}_canonical", con=connection, index=False)
print('Done.')
