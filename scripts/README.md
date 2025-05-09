# CML to LT Converter (cml2lt.py)

This script converts Chemical Markup Language (CML) files to LAMMPS template (LT) files, with special handling for polymer monomers. It ensures consistent atom ordering across related monomer files by using an internal monomer as a reference.

## Requirements

- Python 3.x
- RDKit (installed in conda environment "polymer")
- Required Python packages:
  - rdkit
  - numpy
  - xml.etree.ElementTree

## Installation

1. Ensure you have the conda environment "polymer" activated:
```bash
conda activate polymer
```

2. Verify RDKit is installed in the environment:
```bash
python -c "import rdkit; print(rdkit.__version__)"
```

## Usage

The script requires two main arguments:
- An internal monomer file (must end with 'i.cml')
- A list of monomer files to process

```bash
python cml2lt.py --internal PLAi_T1.cml --monomers PLA_T1.cml PLA_T2.cml
```

Or using the short form:
```bash
python cml2lt.py -i PLAi_T1.cml -m PLA_T1.cml PLA_T2.cml
```

## How It Works

1. **Internal Monomer Processing**:
   - The script first processes the internal monomer file (ending with 'i.cml')
   - Determines the end atom types and their ordering
   - Generates a .lt file for the internal monomer

2. **Other Monomers Processing**:
   - Processes each monomer file in the specified order
   - Uses the end atom types from the internal monomer to ensure consistent ordering
   - Generates .lt files for each monomer

3. **Atom Type Assignment**:
   - Uses RDKit to assign OPLS atom types
   - Ensures consistent atom typing across all monomers

## Output

For each input .cml file, the script generates a corresponding .lt file with:
- Proper atom ordering (left end first, right end second)
- Assigned OPLS atom types
- Bond information
- Coordinates and charges

## Error Handling

The script includes checks for:
- File existence
- Correct internal monomer file naming (*i.cml)
- Valid CML file format
- Successful atom type assignment

## Example

```bash
# Convert a set of monomer files
python cml2lt.py -i PLAi_T1.cml -m PLA_T1.cml PLA_T2.cml

# Output files:
# - PLAi_T1.lt
# - PLA_T1.lt
# - PLA_T2.lt
```

## Notes

- The internal monomer file must end with 'i.cml'
- All input files must be valid CML files
- The script assumes the presence of the OPLS force field definition file
- Generated .lt files are compatible with LAMMPS and moltemplate

## Troubleshooting

If you encounter issues:
1. Verify the conda environment is activated
2. Check that RDKit is properly installed
3. Ensure all input files are valid CML files
4. Verify the internal monomer file ends with 'i.cml'
5. Check file permissions in the working directory
