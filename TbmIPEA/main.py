import pandas as pd
import numpy as np
from contextlib import contextmanager
from utils import *

def getNeutralEnergies():
    with cd('neutral'):
        tot = adiabatic()
        scf = vertical()
    return [scf, tot]


def getIonEnergies():
    ions = {}
    for ion in ['cation', 'anion']:
        with cd(ion):
            with cd('jobex'):
                a = adiabatic()
            with cd('dscf'):
                v = vertical()
            pairs[ion] = {
                    'adiabatic':a,
                    'vertical':v
                    }
    return pairs


def getExcited():
    with cd('ex'):
        with cd('results'):
            with open('exspectrum', 'r') as f:
                lines = f.readlines()
    lines = [line.strip() for line in lines]
    for idx, line in enumerate(lines):
        if '1 a' in line:
            res = line.split()[3]
    return res


def adiabatic():
    with cd('results'):
        with open('job.last', 'r') as f:
            lines = f.readlines()
        lines = [line.strip() for line in lines]
        for idx, line in enumerate(lines):
            if "Total energy + OC corr." in line:
                res = float(line.split()[-1])
    return res


def vertical():
    with cd('results'):
        output = run(['grep', 'cycle', 'gradient'])
        res = output[-1].split()[6]
    return res

def ipea(n, a, c):
    ip = a - n
    ea = n - c
    return ip, ea

def hToEV(h):
    ev = h*27.2114
    return round(ev, 3)

neutrals = getNeutralEnergies()
ionics = getIonEnergies()
ex = getExcited()

for idx, t in enumerate(['vertical', 'adiabatic']):
    n = neutrals[idx]
    a = ionics['anion'][t]
    c = ionics['cation'][t]
    ip, ea = ipea(n, a, c)

    print(f'=> {t.upper()}\n')
    print(f'    IP = {hToEV(ip)}')
    print(f'    EA = {hToEV(ea)}\n')

print(f'=> Optical Gap\n')
print(f'    a1 = {ex}\n')