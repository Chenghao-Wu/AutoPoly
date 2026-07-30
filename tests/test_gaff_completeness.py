"""
GAFF Atom Type Completeness Validation Test

This module validates that all atom types defined in gaff_lt.fdefn
are also present in gaff.lt to ensure proper force field coverage.
"""

import re
import pathlib
from typing import Set, Dict, Tuple
from collections import defaultdict


# File paths relative to project root
GAFF_LT_FDEFN = pathlib.Path(__file__).parent.parent / "AutoPoly/extern/rdlt_data/gaff_lt.fdefn"
GAFF_LT = pathlib.Path(__file__).parent.parent / "AutoPoly/extern/moltemplate/force_fields/gaff.lt"
GAFF_TOMOLTEMPLATE = pathlib.Path(__file__).parent.parent / "AutoPoly/extern/rdlt_data/gaff_tomoltemplate.txt"


def parse_fdefn_atom_types(filepath: pathlib.Path) -> Set[str]:
    """
    Parse gaff_lt.fdefn to extract all @atom:TYPE definitions.

    Args:
        filepath: Path to gaff_lt.fdefn file

    Returns:
        Set of atom type names (e.g., {'c', 'c1', 'c3', ...})
    """
    atom_types = set()
    pattern = re.compile(r'DefineFeature\s+@atom:(\w+)\s+')

    with open(filepath, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.strip().startswith('#'):
                continue
            match = pattern.search(line)
            if match:
                atom_types.add(match.group(1))

    return atom_types


def parse_gaff_lt_atom_types(filepath: pathlib.Path) -> Set[str]:
    """
    Parse gaff.lt to extract all @atom:TYPE declarations from the Data Masses section.

    Args:
        filepath: Path to gaff.lt file

    Returns:
        Set of atom type names (e.g., {'c', 'c1', 'c3', ...})
    """
    atom_types = set()
    pattern = re.compile(r'@atom:(\w+)\s+\d+\.\d+')
    in_masses_section = False

    with open(filepath, 'r') as f:
        for line in f:
            # Look for the Data Masses section
            if 'write_once("Data Masses")' in line:
                in_masses_section = True
                continue
            # End of masses section
            if in_masses_section and line.strip().startswith('}') and 'masses' in line.lower():
                break
            # Extract atom types from masses section
            if in_masses_section:
                match = pattern.search(line)
                if match:
                    atom_types.add(match.group(1))

    return atom_types


def parse_gaff_lt_pair_coeffs(filepath: pathlib.Path) -> Set[str]:
    """
    Parse gaff.lt to extract all @atom:TYPE declarations from pair_coeff section.
    This catches additional atom types that might not be in the masses section.

    Args:
        filepath: Path to gaff.lt file

    Returns:
        Set of atom type names
    """
    atom_types = set()
    pattern = re.compile(r'pair_coeff\s+@atom:(\w+)\s+@atom:')
    in_pair_coeff_section = False

    with open(filepath, 'r') as f:
        for line in f:
            # Look for In Settings section with pair_coeff
            if 'write_once("In Settings")' in line:
                in_pair_coeff_section = True
                continue
            # End of settings section
            if in_pair_coeff_section and line.strip().startswith('}') and 'settings' in line.lower():
                break
            # Extract atom types from pair_coeff lines
            if in_pair_coeff_section and 'pair_coeff' in line:
                # Extract both atom types from the pair_coeff line
                matches = pattern.findall(line)
                for match in matches:
                    atom_types.add(match)

    return atom_types


def get_all_gaff_lt_atom_types(filepath: pathlib.Path) -> Set[str]:
    """
    Get all atom types from gaff.lt by combining masses and pair_coeff sections.

    Args:
        filepath: Path to gaff.lt file

    Returns:
        Set of all atom type names
    """
    masses_types = parse_gaff_lt_atom_types(filepath)
    pair_coeff_types = parse_gaff_lt_pair_coeffs(filepath)
    return masses_types | pair_coeff_types


def parse_mapping_file(filepath: pathlib.Path) -> Dict[str, Dict[str, str]]:
    """
    Parse gaff_tomoltemplate.txt to extract SMARTS patterns and descriptions.

    Args:
        filepath: Path to gaff_tomoltemplate.txt file

    Returns:
        Dictionary mapping atom type symbol to its info:
        {
            'c3': {'smarts': '[CD4]', 'description': 'Sp3 C', 'gaff_name': 'gaff_c3'},
            ...
        }
    """
    mapping = {}
    current_element = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comments and section headers
            if not line or line.startswith('#') or '===' in line:
                continue

            # Parse pipe-delimited lines
            parts = [p.strip() for p in line.split('|')]
            if len(parts) >= 7:
                element, symbol, gaff_name, smarts, moltemplate_type, charge, description = parts[:7]
                # Remove quotes from description
                description = description.strip('"')
                # Clean up smarts pattern
                smarts = smarts.strip()

                if symbol and not symbol.startswith('#'):
                    mapping[symbol] = {
                        'smarts': smarts,
                        'description': description,
                        'gaff_name': gaff_name,
                        'element': element
                    }

    return mapping


def validate_completeness(
    fdefn_types: Set[str],
    gaff_lt_types: Set[str]
) -> Tuple[Set[str], Set[str]]:
    """
    Validate completeness of atom type definitions.

    Args:
        fdefn_types: Atom types from gaff_lt.fdefn
        gaff_lt_types: Atom types from gaff.lt

    Returns:
        Tuple of (missing_in_gaff_lt, extra_in_gaff_lt)
        - missing_in_gaff_lt: Types in fdefn but NOT in gaff.lt (ERROR)
        - extra_in_gaff_lt: Types in gaff.lt but NOT in fdefn (WARNING)
    """
    missing_in_gaff_lt = fdefn_types - gaff_lt_types
    extra_in_gaff_lt = gaff_lt_types - fdefn_types

    return missing_in_gaff_lt, extra_in_gaff_lt


def generate_report(
    fdefn_types: Set[str],
    gaff_lt_types: Set[str],
    missing: Set[str],
    extra: Set[str],
    mapping: Dict[str, Dict[str, str]]
) -> str:
    """
    Generate a structured validation report.

    Args:
        fdefn_types: Atom types from gaff_lt.fdefn
        gaff_lt_types: Atom types from gaff.lt
        missing: Types in fdefn but NOT in gaff.lt
        extra: Types in gaff.lt but NOT in fdefn
        mapping: Mapping file data

    Returns:
        Formatted report string
    """
    lines = [
        "=" * 70,
        "GAFF Atom Type Validation Report",
        "=" * 70,
        "",
        f"Total types in gaff_lt.fdefn: {len(fdefn_types)}",
        f"Total types in gaff.lt: {len(gaff_lt_types)}",
        "",
    ]

    # Report missing types (ERROR)
    if missing:
        lines.extend([
            f"[ERROR] Types in fdefn but MISSING in gaff.lt: {len(missing)}",
            "These types will cause RUNTIME FAILURES because they lack force field parameters:",
            ""
        ])
        for t in sorted(missing):
            desc = mapping.get(t, {}).get('description', 'No description')
            lines.append(f"  - @atom:{t} ({desc})")
        lines.append("")
    else:
        lines.extend([
            "[OK] All atom types in gaff_lt.fdefn exist in gaff.lt",
            ""
        ])

    # Report extra types (WARNING)
    if extra:
        lines.extend([
            f"[WARNING] Types in gaff.lt but NOT in fdefn: {len(extra)}",
            "These types are valid in GAFF but unused in current fdefn:",
            ""
        ])
        for t in sorted(extra):
            desc = mapping.get(t, {}).get('description', 'No description')
            lines.append(f"  - @atom:{t} ({desc})")
        lines.append("")
    else:
        lines.extend([
            "[OK] All atom types in gaff.lt are defined in fdefn",
            ""
        ])

    # Final verdict
    lines.append("-" * 70)
    if missing:
        lines.append("Validation: FAILED - Missing force field parameters")
    else:
        lines.append("Validation: PASSED - All types have force field parameters")
    lines.append("=" * 70)

    return "\n".join(lines)


def test_gaff_completeness():
    """
    Main test function for pytest.

    This test validates that all atom types in gaff_lt.fdefn have
    corresponding definitions in gaff.lt.
    """
    # Parse files
    fdefn_types = parse_fdefn_atom_types(GAFF_LT_FDEFN)
    gaff_lt_types = get_all_gaff_lt_atom_types(GAFF_LT)
    mapping = parse_mapping_file(GAFF_TOMOLTEMPLATE)

    # Validate completeness
    missing, extra = validate_completeness(fdefn_types, gaff_lt_types)

    # Generate and print report
    report = generate_report(fdefn_types, gaff_lt_types, missing, extra, mapping)
    print(report)

    # Assert that there are no missing types
    assert not missing, (
        f"Found {len(missing)} atom types in gaff_lt.fdefn that are missing in gaff.lt. "
        "These will cause runtime failures. See report above for details."
    )


def test_fdefn_types_exist():
    """Basic smoke test to ensure fdefn file is being parsed correctly."""
    fdefn_types = parse_fdefn_atom_types(GAFF_LT_FDEFN)
    # We should have at least some common types
    assert 'c3' in fdefn_types, "Expected @atom:c3 to be in fdefn"
    assert 'ha' in fdefn_types, "Expected @atom:ha to be in fdefn"
    assert 'oh' in fdefn_types, "Expected @atom:oh to be in fdefn"
    assert len(fdefn_types) > 20, "Expected more atom types in fdefn"


def test_gaff_lt_types_exist():
    """Basic smoke test to ensure gaff.lt file is being parsed correctly."""
    gaff_lt_types = get_all_gaff_lt_atom_types(GAFF_LT)
    # We should have many more types in gaff.lt
    assert 'c3' in gaff_lt_types, "Expected @atom:c3 to be in gaff.lt"
    assert 'ha' in gaff_lt_types, "Expected @atom:ha to be in gaff.lt"
    assert len(gaff_lt_types) > 60, "Expected more atom types in gaff.lt"


if __name__ == "__main__":
    # Run validation when executed directly
    print(f"Validating GAFF atom type definitions...\n")
    print(f"gaff_lt.fdefn: {GAFF_LT_FDEFN}")
    print(f"gaff.lt: {GAFF_LT}")
    print(f"gaff_tomoltemplate.txt: {GAFF_TOMOLTEMPLATE}\n")

    test_gaff_completeness()
