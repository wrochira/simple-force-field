"""
Basic structural tests for force field definitions.

This test module performs fundamental sanity checks on force field classes,
verifying that:
- All force fields define atom types and interaction types
- Parameters and constants match the potential function definitions
- Force (F) and potential energy (U) values are finite and sensible
- Cutoff behaviour is correctly implemented (F and U vanish beyond cutoff)
- Potential energy is continuous at the cutoff distance

The tests are designed to catch structural issues and basic implementation
errors without requiring detailed knowledge of the underlying physics.
"""

import sys
from pathlib import Path
from math import isclose, isfinite

# Add parent directory to path to import force_field module
sys.path.insert(0, str(Path(__file__).parent.parent))

from force_field import (
    FBA,
    OPC3,
    SPC,
    SPC_E,
    SPC_FW,
    TIP3P,
    TIP3P_CHARMM,
    TIP3P_EW,
    TIP3P_FB,
    TIP3P_FS,
    TIP3P_PW,
)

FORCE_FIELDS = [
    FBA,
    OPC3,
    SPC,
    SPC_E,
    SPC_FW,
    TIP3P,
    TIP3P_CHARMM,
    TIP3P_EW,
    TIP3P_FB,
    TIP3P_FS,
    TIP3P_PW,
]


def check_interaction_basic(interaction):
    """Basic structural checks on an interaction type."""
    # Parameter and constant dictionaries should match the potential definition
    potential = interaction.potential
    assert set(interaction.parameters.keys()) == set(
        potential.parameter_names
    ), f'Parameters for {interaction.name} do not match potential {potential.name}'
    assert set(interaction.constants.keys()) == set(
        potential.constant_names
    ), f'Constants for {interaction.name} do not match potential {potential.name}'

    # Must involve at least two atom types
    assert len(interaction.atom_types) >= 2, f'{interaction.name} has too few atom types'


def check_interaction_values(interaction):
    """Check that F and U behave sensibly around the cutoff."""
    cutoff = interaction.cutoff

    # Choose a few representative distances safely away from singularities
    test_distances = [0.8, 0.9, 1.0, 1.1, 1.2]

    for r in test_distances:
        F = interaction.F(r)
        U = interaction.U(r)

        # Basic sanity: values should be finite
        assert isfinite(F), f'F({r}) is not finite for {interaction.name}'
        assert isfinite(U), f'U({r}) is not finite for {interaction.name}'

        # Cutoff behaviour
        if cutoff is not None and r > cutoff:
            assert F == 0.0, f'F({r}) != 0 beyond cutoff for {interaction.name}'
            assert U == 0.0, f'U({r}) != 0 beyond cutoff for {interaction.name}'

    # At the cutoff, U should be continuous and zero by construction
    if cutoff is not None:
        U_cut = interaction.U(cutoff)
        assert isclose(U_cut, 0.0, abs_tol=1e-8), (
            f'U(cutoff) != 0 for {interaction.name}: U({cutoff}) = {U_cut}'
        )


def test_force_field_class(ff_type):
    """Run a small battery of checks on a force-field class."""
    print(f'\n=== Testing {ff_type.__name__} ===')
    ff = ff_type()

    # Force field should expose at least one atom and one interaction type
    assert ff.atom_types, f'{ff_type.__name__} defines no atom types'
    assert ff.interaction_types, f'{ff_type.__name__} defines no interaction types'

    for interaction in ff.interaction_types:
        print(f'  - {interaction.name}')
        check_interaction_basic(interaction)
        check_interaction_values(interaction)

    print(f'{ff_type.__name__}: OK')


if __name__ == '__main__':
    for ff_type in FORCE_FIELDS:
        test_force_field_class(ff_type)
