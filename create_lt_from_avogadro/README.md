# CML to LT File Converter for AutoPoly

This module converts CML (Chemical Markup Language) files to LT (LAMMPS Template) files for use with the AutoPoly polymerization system.

## Overview

The converter generates three types of monomer .lt files required for polymer building:
- **Internal monomer** (suffix 'i'): Two connection points for linking to other monomers
- **Left-end monomer** (suffix 'le'): One terminal end and one connection point  
- **Right-end monomer** (suffix 're'): One connection point and one terminal end

## Features

- ✅ **CML File Parsing**: Parses CML files with atom coordinates and bond information
- ✅ **Atom Type Assignment**: Assigns correct OPLS atom types for polymer connectivity
- ✅ **3D Coordinate Preservation**: Maintains exact 3D coordinates from CML files
- ✅ **Bond Connectivity**: Preserves all bond information for proper molecular structure
- ✅ **AutoPoly Compatibility**: Generates .lt files compatible with the AutoPoly system
- ✅ **LOPLS Support**: Optional LOPLS force field support
- ✅ **Error Handling**: Robust error handling for missing files and parsing errors

## Atom Type Mapping

| Monomer Type | C1 Atom Type | C2 Atom Type | Description |
|--------------|--------------|--------------|-------------|
| Internal (i) | @atom:81 | @atom:82 | Both connection points (CH2) |
| Left-end (le) | @atom:80 | @atom:82 | Terminal (CH3) + connection (CH2) |
| Right-end (re) | @atom:81 | @atom:81 | Connection (CH2) + terminal (CH2) |

## Installation

No additional installation required beyond the existing AutoPoly dependencies:
- Python 3.x
- RDKit
- XML parsing (built-in)

## Usage

### Command Line Interface

```bash
# Basic usage
python rdlt_avogadro.py --internal PEi.cml --left PEl.cml --right PEr.cml --base-name PE --output-dir output

# With LOPLS force field
python rdlt_avogadro.py --internal PEi.cml --left PEl.cml --right PEr.cml --base-name PE --output-dir output --lopls

# With custom feature definitions
python rdlt_avogadro.py --internal PEi.cml --left PEl.cml --right PEr.cml --base-name PE --output-dir output --opls-fdef path/to/opls_lt.fdefn
```

### Python API

```python
from rdlt_avogadro import CMLToLTConverter

# Create converter
converter = CMLToLTConverter()

# Process all three monomer variants
cml_files = {
    'internal': 'PEi.cml',
    'left': 'PEl.cml', 
    'right': 'PEr.cml'
}

generated_files = converter.process_polymer_monomers(cml_files, 'PE', 'output')

# Process single monomer
converter.process_cml_to_lt('PEi.cml', 'PEi', 'internal', 'PEi.lt')
```

## Testing

The implementation has been thoroughly tested with the provided PE (Polyethylene) CML files:

### Test Results

✅ **Command Line Interface**: Works correctly with all options
✅ **File Generation**: All three monomer variants generated successfully
✅ **Atom Type Assignment**: Correct OPLS atom types assigned for each variant
✅ **Coordinate Preservation**: 3D coordinates preserved from CML files
✅ **Bond Connectivity**: All bonds maintained correctly
✅ **File Structure**: Generated files have correct AutoPoly-compatible structure
✅ **Error Handling**: Gracefully handles missing files and parsing errors
✅ **LOPLS Support**: LOPLS flag adds appropriate import statements

### Generated Files

For PE monomers, the following files are generated:
- `PEi.lt`: Internal monomer with 6 atoms (2 C, 4 H) and 5 bonds
- `PEle.lt`: Left-end monomer with 7 atoms (2 C, 5 H) and 6 bonds  
- `PEre.lt`: Right-end monomer with 7 atoms (2 C, 5 H) and 6 bonds

### File Structure

Each generated .lt file contains:
```lt
import "oplsaa.lt"    # <-- defines the standard "OPLSAA" force field
{monomer_name} inherits OPLSAA {

# atom-id  mol-id  atom-type charge      X         Y        Z
  write("Data Atoms") {
    $atom:C1 $mol:... @atom:81 0.00    0.136    0.013   -0.007
    $atom:C2 $mol:... @atom:82 0.00    1.652   -0.002    0.014
    ...
  }

  write('Data Bond List') {
    $bond:C1C2	$atom:C1	$atom:C2
    ...
  }
} # {monomer_name}
```

## Integration with AutoPoly

The generated .lt files are fully compatible with the AutoPoly polymerization system:

```python
import AutoPoly

# Use the generated monomers
system = AutoPoly.System(out="test_polymer")
polymer = AutoPoly.Polymer(ChainNum=10, Sequence=["PEi"]*50, topology="linear")
poly = AutoPoly.Polymerization(name="TestPE", system=system, model=[polymer], run=True)
```

## Example Files

The repository includes example CML files for testing:
- `PEi.cml`: Internal polyethylene monomer
- `PEl.cml`: Left-end polyethylene monomer  
- `PEr.cml`: Right-end polyethylene monomer

## Test Scripts

- `test_cml_to_lt.py`: Basic functionality tests
- `example_usage.py`: Usage examples
- `integration_test.py`: AutoPoly integration tests

## Error Handling

The converter handles various error conditions:
- Missing CML files: Skips with warning
- Invalid CML format: Raises descriptive error
- Missing atom types: Validates before file generation
- File write errors: Handles gracefully

## Performance

- Fast parsing of CML files
- Efficient RDKit molecule building
- Minimal memory usage
- Quick .lt file generation

## Future Enhancements

Potential improvements:
- Support for more complex molecular structures
- Automatic connection point detection
- Integration with molecular visualization tools
- Batch processing of multiple monomer sets
- Support for additional force fields

## License

This module is part of the AutoPoly project and follows the same licensing terms. 