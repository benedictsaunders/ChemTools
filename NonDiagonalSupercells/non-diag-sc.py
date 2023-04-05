import numpy as np
from sys import argv
from numpy.linalg import inv, det
from ase.build import make_supercell, niggli_reduce, bulk
from pprint import pprint
from functools import reduce
from math import ceil

# Given n, return a list of lists which contains 3 integers that multiply to produce that number
def get_triple_factors(n):
    triples = []
    def factors(x):
        pairs = []
        fs = list(set(reduce(list.__add__, ([i, x//i] for i in range(1, int(x**0.5) + 1) if x % i == 0))))
        fs.sort()
        for i in range(0, ceil(len(fs)/2)):
            pairs.append(
                    [
                        fs[i],
                        fs[-1-i]
                    ]
            )
        return pairs
    for pair in factors(n):
        for p in pair:
            pfactors = factors(p)




def non_diagonal(atoms, n):
    def conditions(S):
        return (0 <= S[0][1] and S[1][1] and 0 <= S[0][2] and S[1][2] < S[2][2])
    #niggli_reduce
    niggli_reduce(atoms)
    print(atoms)

    S = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 2]
    ])

    # make sure that we formed the utHNF correctly
    assert(int(S[0][0]*S[1][1]*S[2][2]) == int(np.linalg.det(S)))


    print(det(S))
    #Rules:
    
        # 0 <= S12 < S22
        # 0 < S13 , S23 < S33
    
        # Adding integer multiple of one row to another row
        # Interchanging two rows
        # Multiplying a row by -1

    # prod = S11*S22*S33 = detS
    if conditions(S):
        return make_supercell(atoms, S)

def makeNDSM(n):
    """
    Based on the algorithm by Lloyd-Williams and Monserratt in doi:10.1103/PhysRevB.92.184301
        Parameters:
            n int: Number of primitives to be placed in the superecell.
    """

    count_hnf = 0
    for a in range(1, n+1):
        if not (n%a == 0):
            continue
        quotient = int(n/a)
        for c in range(1, quotient+1):
            if not (quotient%c == 0):
                continue
            f = quotient/c
            count_hnf=count_hnf+c*f**2
            print(count_hnf)
    num_hnf = count_hnf
    count_hnf = 0

    hnf = np.zeros(
        (int(num_hnf),3,3)
    )

    for a in range(1, n+1):
        if not (n%a == 0):
            continue
        quotient = int(n/a)
        for c in range(1, quotient+1):
            if not (quotient%c == 0):
                continue
            f = int(quotient/c)
            for b in range(0, c):
                for d in range(0, f):
                    for e in range(0, f):
                        hnf[count_hnf][0][0] = a
                        hnf[count_hnf][0][1] = b
                        hnf[count_hnf][1][1] = c
                        hnf[count_hnf][0][2] = d
                        hnf[count_hnf][1][2] = e
                        hnf[count_hnf][2][2] = f
                        count_hnf += 1

    return hnf
if __name__ == "__main__":
   m =  makeNDSM(int(argv[1]))
   pprint(m)