#!/usr/bin/env python3
"""
Example: Programmatic usage of CMLToLTConverterTacticity with PDMS monomers
"""
from rdlt_avogadro import CMLToLTConverterTacticity

# Optionally, define a custom atom type mapping (edit as needed)
custom_mapping = {
    'internal': {'O1': '@atom:122', 'Si2': '@atom:80'},
    'left': {'O1': '@atom:122', 'Si2': '@atom:80'},
    'right': {'O1': '@atom:122', 'Si2': '@atom:80'}
}

# Create the converter (with or without custom mapping)
converter = CMLToLTConverterTacticity(atom_type_mapping=custom_mapping)

# Define your CML files for PDMS (make sure these files exist in the current directory)
cml_files = {
    'internal': 'PDMSi.cml',
    'left': 'PDMSle.cml',
    'right': 'PDMSre.cml'
}

# Generate tacticity variants (e.g., atactic)
generated_files = converter.generate_tacticity_variants(
    cml_files, base_name='PDMS', output_dir='PDMS'
)

print("Generated files:")
for monomer_type, files in generated_files.items():
    print(f"{monomer_type}: {files}") 