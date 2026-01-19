#!/usr/bin/env python3
"""
Unit tests for Gasteiger charge calculation in LTWriter.

Tests that Gasteiger charges are automatically calculated for GAFF force field.
"""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from AutoPoly.monomer_generator import LTWriter, SMARTSTyper, MonomerVariant


def test_gasteiger_charges_calculated():
    """Test that Gasteiger charges are calculated for GAFF."""
    # Create simple molecule
    mol = Chem.MolFromSmiles('CC')
    mol = Chem.AddHs(mol)

    # Generate 3D conformer
    AllChem.EmbedMolecule(mol, randomSeed=42)

    # Assign atom types
    typer = SMARTSTyper('gaff')
    typer.assign_atom_types(mol)

    # Create LTWriter with empty charge dict (like GAFF)
    writer = LTWriter('gaff', {}, verbose=True)

    # Ensure charges are calculated
    writer._ensure_gasteiger_charges(mol)

    # Verify charges were calculated and cached
    assert len(writer.charge_dict) > 0, "No charges were calculated"

    # Verify at least some non-zero charges
    non_zero_charges = [charge for charge in writer.charge_dict.values() if charge != 0.0]
    assert len(non_zero_charges) > 0, "All charges are zero"


def test_gasteiger_neutral_molecule():
    """Test that Gasteiger charges sum to ~0.0 for neutral molecule."""
    mol = Chem.MolFromSmiles('CC')
    mol = Chem.AddHs(mol)

    # Generate 3D conformer
    AllChem.EmbedMolecule(mol, randomSeed=42)

    typer = SMARTSTyper('gaff')
    typer.assign_atom_types(mol)

    writer = LTWriter('gaff', {})
    writer._ensure_gasteiger_charges(mol)

    # Calculate net charge
    total_charge = sum(writer.charge_dict.values())

    # Net charge should be close to 0.0 (Gasteiger is approximate)
    assert abs(total_charge) < 0.05, f"Net charge {total_charge} is not close to 0.0"


def test_opls_unchanged():
    """Test that OPLS force field is unchanged."""
    mol = Chem.MolFromSmiles('CC')
    mol = Chem.AddHs(mol)

    # Generate 3D conformer
    AllChem.EmbedMolecule(mol, randomSeed=42)

    typer = SMARTSTyper('oplsaa')
    typer.assign_atom_types(mol)

    # Pre-populate charge dict (like OPLS database)
    initial_charges = {'@atom:80': -0.18, '@atom:81': 0.06}
    writer = LTWriter('oplsaa', initial_charges.copy())

    writer._ensure_gasteiger_charges(mol)

    # Should not add new entries (OPLS uses database)
    assert len(writer.charge_dict) == 2, "OPLS should not calculate Gasteiger charges"
    assert writer.charge_dict == initial_charges, "OPLS charge dict should be unchanged"


def test_gasteiger_caching():
    """Test that charges are cached and not recalculated."""
    mol = Chem.MolFromSmiles('CC')
    mol = Chem.AddHs(mol)

    # Generate 3D conformer
    AllChem.EmbedMolecule(mol, randomSeed=42)

    typer = SMARTSTyper('gaff')
    typer.assign_atom_types(mol)

    writer = LTWriter('gaff', {})

    # First call - should calculate
    writer._ensure_gasteiger_charges(mol)
    first_dict = writer.charge_dict.copy()

    # Second call - should use cached values
    writer._ensure_gasteiger_charges(mol)
    second_dict = writer.charge_dict.copy()

    assert first_dict == second_dict, "Charges should be cached"


def test_gasteiger_invalid_values():
    """Test that invalid Gasteiger charges are handled gracefully."""
    mol = Chem.MolFromSmiles('C[SiH2]C')  # Silane - can produce unusual charges
    mol = Chem.AddHs(mol)

    # Generate 3D conformer
    AllChem.EmbedMolecule(mol, randomSeed=42)

    typer = SMARTSTyper('gaff')
    typer.assign_atom_types(mol)

    writer = LTWriter('gaff', {})

    # Should not raise exception
    try:
        writer._ensure_gasteiger_charges(mol)
        # Verify all charges are valid (not NaN or Inf)
        for atom_type, charge in writer.charge_dict.items():
            assert charge == charge, f"NaN charge found for {atom_type}"
            assert abs(charge) < 1000, f"Suspiciously large charge for {atom_type}: {charge}"
    except Exception as e:
        pytest.fail(f"Gasteiger calculation raised exception: {e}")


def test_gasteiger_pre_existing_charges():
    """Test that pre-existing charges are preserved."""
    mol = Chem.MolFromSmiles('CC')
    mol = Chem.AddHs(mol)

    # Generate 3D conformer
    AllChem.EmbedMolecule(mol, randomSeed=42)

    typer = SMARTSTyper('gaff')
    typer.assign_atom_types(mol)

    # Pre-populate with some charges
    existing_charge = -0.5
    writer = LTWriter('gaff', {'@atom:ct': existing_charge})

    writer._ensure_gasteiger_charges(mol)

    # Existing charge should be preserved
    assert writer.charge_dict.get('@atom:ct') == existing_charge, \
        "Pre-existing charges should not be overwritten"


def test_write_variant_with_gasteiger():
    """Integration test: write variant with Gasteiger charges."""
    mol = Chem.MolFromSmiles('CC')  # Ethane - simple molecule without dummy atoms
    mol = Chem.AddHs(mol)

    # Generate 3D conformer
    AllChem.EmbedMolecule(mol, randomSeed=42)

    typer = SMARTSTyper('gaff')
    typer.assign_atom_types(mol)

    # Create variant
    variant = MonomerVariant(
        base_name="test",
        variant_type="single",
        mol=mol,
        smiles=Chem.MolToSmiles(mol),
        connection_atoms=(0, 1),  # Both carbons
        force_field="gaff",
        position=0,
        atom_ids={}
    )

    writer = LTWriter('gaff', {})

    # Write to string
    import io
    output = io.StringIO()

    # This should trigger Gasteiger calculation
    writer._write_atoms_block(output, variant)

    output_str = output.getvalue()

    # Verify output contains atoms
    assert '$atom:' in output_str, "Output should contain atom definitions"

    # Verify charges were calculated
    assert len(writer.charge_dict) > 0, "Charges should be calculated"

    # Verify not all charges are zero (should work for ethane)
    non_zero = [c for c in writer.charge_dict.values() if c != 0.0]
    assert len(non_zero) > 0, "Some charges should be non-zero for ethane"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
