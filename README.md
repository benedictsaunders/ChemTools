# ConformerSearch
This tool is able to generate conformers from a given smile string. All of the conformers are then minimised with MMFF94, and the 5 geometries with the lowest forcefield energy are outputted as XYZ files in the form `[min_energy_order]_[conformer_ID].xyz`. However, all the geomteries found, both relaxed and not, can be accessed too. Before and after MMFF94 optmisation, the RDKit `mol` objects into which the conformers are embedded are [pickled](https://docs.python.org/3/library/pickle.html) and saved as `.pickle` files.

## Requirements
This tools requires the following python libraries:
* RDkit
* Pickle
* Argparse
* MatPlotLib
* Pandas

## Usage 

The tool can be launched from the command line, `python -s <smiles> -n <number of conformers>` for a single molecule.

For multiple molecule, one can run the supplied bash script. The list of SMILES can be given in a file `smiles-list.txt`. Each molecule will be saved in its own directory, named numerically. To name the directories differently, you can provide a second file, `names-list.txt`, a list of the desired names of the desired directories.
