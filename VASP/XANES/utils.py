from contextlib import contextmanager
import os
from pprint import pprint

class AliasDict(dict):
    """
    A derivative of the dict class to allow for dictionary keys to be aliased.
    Stolen from jasonharper on StackOverflow (I won't tel; if you don't!)
    Thanks Jason :)
    """
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self.aliases = {}

    def __getitem__(self, key):
        return dict.__getitem__(self, self.aliases.get(key, key))

    def __setitem__(self, key, value):
        return dict.__setitem__(self, self.aliases.get(key, key), value)

    def add_alias(self, key, alias):
        self.aliases[alias] = key

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    if not os.path.exists(f"{prevdir}/{newdir}"):
        os.mkdir(newdir)
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def handle_magmoms(magmoms, symbols):
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
    elements = magmoms[::2]
    values = magmoms[1::2]
    d = dict(zip(elements, values))
    for s in symbols:
        if s not in d.keys():
            d[s] = 0
    return d

def handle_hubbard(luj, symbols):
    """
    Takes the argparse input and transform it to a nested dictionary, 
    and accounts for unmentioed species in the correction.
    """
    if luj is None:
        print("Hubbard corrections not set.")
        return None
    labels = ["l", "U", "J"]
    n = 4
    elements = []
    d = AliasDict()
    separated = [luj[i: i + n] for i in range(0, len(luj), n)]
    for indiv in separated:
        elements.append(indiv[0])
        d[indiv[0]] = AliasDict(zip(labels, [float(x) for x in indiv[-3:]]))
    for s in symbols:
        if s not in elements:
            d[s] = AliasDict(
                zip(labels, [-1, 0, 0])
            )  # Negative 1 for no onsite interation added.
    return d
