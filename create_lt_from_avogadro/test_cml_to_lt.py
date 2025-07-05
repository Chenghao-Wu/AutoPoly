#!/usr/bin/env python3
"""
Test script for CML to LT conversion.
"""

import os
import sys
from pathlib import Path

# Add the current directory to the path
sys.path.insert(0, str(Path(__file__).parent))

from rdlt_avogadro import CMLToLTConverter

def test_pe_monomers():
    """Test conversion of PE (Polyethylene) CML files."""
    
    print("Testing CML to LT conversion for PE monomers...")
    
    # Define CML files
    cml_files = {
        'internal': 'PEi.cml',
        'left': 'PEl.cml',
        'right': 'PEr.cml'
    }
    
    # Check if files exist
    for monomer_type, cml_file in cml_files.items():
        if not os.path.exists(cml_file):
            print(f"Error: {cml_file} not found!")
            return
    
    # Create converter
    converter = CMLToLTConverter()
    
    # Process files
    try:
        generated_files = converter.process_polymer_monomers(
            cml_files, 'PE', 'output', loplsflag=False
        )
        
        print("\nSuccess! Generated files:")
        for monomer_type, file_path in generated_files.items():
            if file_path and os.path.exists(file_path):
                print(f"  {monomer_type}: {file_path}")
                
                # Check file content
                with open(file_path, 'r') as f:
                    content = f.read()
                    print(f"    File size: {len(content)} characters")
                    print(f"    Contains 'inherits OPLSAA': {'inherits OPLSAA' in content}")
                    print(f"    Contains 'Data Atoms': {'Data Atoms' in content}")
                    print(f"    Contains 'Data Bond List': {'Data Bond List' in content}")
                    
                    # Check atom types
                    if monomer_type == 'internal':
                        if '@atom:81' in content and '@atom:82' in content:
                            print(f"    ✓ Correct atom types for internal monomer")
                        else:
                            print(f"    ✗ Incorrect atom types for internal monomer")
                    elif monomer_type == 'left':
                        if '@atom:80' in content and '@atom:82' in content:
                            print(f"    ✓ Correct atom types for left-end monomer")
                        else:
                            print(f"    ✗ Incorrect atom types for left-end monomer")
                    elif monomer_type == 'right':
                        if '@atom:81' in content and '@atom:81' in content:
                            print(f"    ✓ Correct atom types for right-end monomer")
                        else:
                            print(f"    ✗ Incorrect atom types for right-end monomer")
            else:
                print(f"  {monomer_type}: FAILED")
                
    except Exception as e:
        print(f"Error: {e}")

def test_single_monomer():
    """Test conversion of a single CML file."""
    
    print("\nTesting single monomer conversion...")
    
    if not os.path.exists('PEi.cml'):
        print("Error: PEi.cml not found!")
        return
    
    converter = CMLToLTConverter()
    
    try:
        output_file = converter.process_cml_to_lt(
            'PEi.cml', 'PEi_test', 'internal', 'PEi_test.lt'
        )
        print(f"Generated single monomer: {output_file}")
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    test_pe_monomers()
    test_single_monomer() 