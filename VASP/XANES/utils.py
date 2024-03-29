from contextlib import contextmanager
import os
from pprint import pprint
import subprocess as sp
from io import TextIOWrapper


class AliasDict(dict):
    """
    A derivative of the dict class to allow for dictionary keys to be aliased.
    Stolen from jasonharper on StackOverflow (I won't tell if you don't!)
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


def get_lines(file_object):
    assert isinstance(file_object, TextIOWrapper)
    return [l.strip() for l in file_object.readlines()]


class input_file:
    def __init__(self, fname) -> None:
        with open(fname, "r") as f:
            for line in get_lines(f):
                pass
        pass


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
    types = [int, float, float]
    n = 4
    elements = []
    d = AliasDict()
    separated = [luj[i : i + n] for i in range(0, len(luj), n)]
    for indiv in separated:
        elements.append(indiv[0])
        d[indiv[0]] = AliasDict(
            zip(labels, [types[i](x) for i, x in enumerate(indiv[-3:])])
        )
    for s in symbols:
        if s not in elements:
            d[s] = AliasDict(
                zip(labels, [-1, 0, 0])
            )  # Negative 1 for no onsite interation added.
    return d


def newline(s):
    if type(s) == list:
        l = [str(x) for x in s]
        s = " ".join(l)
    return str(s) + "\n"


def runcmd(cmd, outlocation=sp.PIPE, errlocation=sp.PIPE):
    if not isinstance(cmd, list):
        cmd = cmd.split()
    output = sp.run(cmd, stdout=outlocation, stderr=errlocation)
    return output
