Script to automate assignment of atom types and build `.lt` files for use with `moltemplate` and the OPLS and L-OPLS forcefields that come packaged with it. The script uses RDKit and its feature definition format. Style of files generated is intended to mimic what was being built by hand.

To script requires a valid SMILES string as an input and prints the generated template file.

To build a template file for hexane:

`python rdlt.py --smi 'CCCCCC' -n HEX -l -c > hex.lt`

The `-l` flag tells the script to use the L-OPLS forcefield to overwrite the standard OPLSAA atom type assignments. The `-c` flag prints the expected atom charges at the end, and sums them. A warning is then printed if they do not sum approximately to zero.

Use the `-h` flag for help with the other options.
