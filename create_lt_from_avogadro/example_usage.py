#!/usr/bin/env python3
"""
Example usage of CML to LT converter.
"""

from rdlt_avogadro import CMLToLTConverter

def main():
    """Example usage of the CML to LT converter."""
    
    # Create converter instance
    converter = CMLToLTConverter()
    
    # Example 1: Process PE (Polyethylene) monomers
    print("Example 1: Processing PE monomers...")
    
    pe_cml_files = {
        'internal': 'PEi.cml',
        'left': 'PEl.cml',
        'right': 'PEr.cml'
    }
    
    try:
        generated_files = converter.process_polymer_monomers(
            pe_cml_files, 'PE', 'output_monomers'
        )
        
        print("Generated PE monomer files:")
        for monomer_type, file_path in generated_files.items():
            if file_path:
                print(f"  {monomer_type}: {file_path}")
        
    except Exception as e:
        print(f"Error processing PE monomers: {e}")
    
    # Example 2: Process single monomer
    print("\nExample 2: Processing single monomer...")
    
    try:
        output_file = converter.process_cml_to_lt(
            'PEi.cml', 'PEi_single', 'internal', 'PEi_single.lt'
        )
        print(f"Generated single monomer: {output_file}")
        
    except Exception as e:
        print(f"Error processing single monomer: {e}")
    
    # Example 3: Using with custom feature definitions
    print("\nExample 3: Using custom feature definitions...")
    
    # You can specify custom feature definition files if available
    # converter = CMLToLTConverter(
    #     opls_fdef_path='path/to/opls_lt.fdefn',
    #     lopls_fdef_path='path/to/lopls_lt.fdefn'
    # )

if __name__ == "__main__":
    main() 