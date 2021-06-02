from contextlib import contextmanager
import os
import subprocess as sp
from time import sleep

def run(command):
    sleep(1)
    p = sp.Popen(command, stdout=sp.PIPE, encoding="utf-8")
    output, errors = p.communicate()
    return output.splitlines()

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
