"""
Comprehensive physical validation tests for interaction potentials.

This test module performs detailed physics-based validation of force field
interaction potentials, including:

- **Harmonic potentials** (bonds/angles):
  * Force and energy vanish at equilibrium
  * Symmetry around equilibrium position
  * Correct restoring force direction
  * Water-specific validation (O-H bond ~0.1 nm, H-O-H angle ~104.5°)

- **Coulombic potentials**:
  * Correct sign based on charge product
  * Magnitude scales correctly with distance (1/r behaviour)
  * Energy increases dramatically at short distances

- **Lennard-Jones potentials**:
  * Minimum at r = 2^(1/6) * sigma
  * Repulsive at short range (U > 0 for r < sigma)
  * Attractive beyond minimum (U < 0)
  * Energy approaches zero at large distances

- **General validation**:
  * Smooth decay towards cutoff distance
  * Exact zeros beyond cutoff
  * All values remain finite

The tests sample each interaction over physically relevant distance ranges
and validate that the potential functions exhibit expected behaviour for
molecular dynamics simulations of water molecules.
"""

import sys
from pathlib import Path
from math import isfinite, isclose

# Add parent directory to path to import force_field module
sys.path.insert(0, str(Path(__file__).parent.parent))

from force_field.potentials import Harmonic, Coulombic, LennardJones

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


def get_distance_range(interaction, default_min: float = 0.1, default_max: float = 1.0) -> tuple[float, float]:
    """Heuristically choose a sensible distance/angle range for an interaction."""
    params = interaction.parameters
    cutoff = interaction.cutoff

    r_min = default_min
    r_max = default_max

    # Prefer explicit equilibrium distance if available (e.g. harmonic bonds / angles)
    if 'eq' in params:
        eq = params['eq']
        r_min = 0.6 * eq
        r_max = 1.4 * eq

    # For Lennard-Jones style interactions, use sigma if available
    elif 'sigma' in params:
        sigma = params['sigma']
        r_min = 0.6 * sigma
        r_max = 2.5 * sigma

    # Otherwise, if a cutoff is defined, use a fraction of it
    elif cutoff is not None:
        r_min = 0.2 * cutoff
        r_max = 1.2 * cutoff

    # Ensure strictly positive and sane ordering
    r_min = max(r_min, 1e-3)
    if r_max <= r_min:
        r_max = r_min * 10.0

    return r_min, r_max


def sample_interaction(interaction, num_points: int = 200):
    """Sample F(x) and U(x) over a range for a single interaction."""
    x_min, x_max = get_distance_range(interaction)
    xs = [x_min + (x_max - x_min) * i / (num_points - 1) for i in range(num_points)]

    data = []
    for x in xs:
        F_val = interaction.F(x)
        U_val = interaction.U(x)
        data.append((x, F_val, U_val))

    return data


def _symbols(interaction):
    return [atom_type.symbol for atom_type in interaction.atom_types]


def is_oh_bond(interaction) -> bool:
    """Detect O-H bond interactions in water models."""
    if interaction.category_name != 'Bond Stretching':
        return False
    syms = sorted(_symbols(interaction))
    return syms == ['HW', 'OW']


def is_hoh_angle(interaction) -> bool:
    """Detect H-O-H angle interactions in water models."""
    if interaction.category_name != 'Angle Bending':
        return False
    syms = _symbols(interaction)
    return len(syms) == 3 and syms[0] == 'HW' and syms[1] == 'OW' and syms[2] == 'HW'


def sense_check_common(interaction, data):
    """Checks that should hold for any interaction."""
    cutoff = interaction.cutoff

    # Finite values across the sampled range
    for x, F_val, U_val in data:
        assert isfinite(F_val), f'{interaction.name}: F({x}) is not finite'
        assert isfinite(U_val), f'{interaction.name}: U({x}) is not finite'

    # Cutoff behaviour:
    # - Beyond cutoff: U,F must be ~0 (this is enforced exactly in InteractionType).
    # - Around cutoff: U should be smoothly tending towards 0; we only enforce this
    #   loosely because we may not sample the cutoff point exactly.
    if cutoff is not None:
        # Check monotonic decay in |U| as we approach the cutoff from below:
        below = sorted([d for d in data if d[0] < cutoff], key=lambda d: d[0])
        if len(below) >= 3:
            # Take three closest points below the cutoff
            tail = below[-3:]
            mags = [abs(d[2]) for d in tail]
            # We expect the magnitude to be non-increasing as we approach the cutoff
            assert mags[0] >= mags[1] >= mags[2] or isclose(mags[0], mags[2], rel_tol=1e-1), (
                f'{interaction.name}: |U| does not decay smoothly towards cutoff; '
                f'|U| near cutoff = {mags}'
            )

        beyond = [d for d in data if d[0] >= cutoff * 1.01]
        for x, F_val, U_val in beyond:
            assert isclose(F_val, 0.0, abs_tol=1e-9), (
                f'{interaction.name}: F(x={x:.3f}) beyond cutoff not 0, got {F_val}'
            )
            assert isclose(U_val, 0.0, abs_tol=1e-9), (
                f'{interaction.name}: U(x={x:.3f}) beyond cutoff not 0, got {U_val}'
            )


def sense_check_harmonic(interaction):
    """Sense-check for harmonic bonds/angles, with water-specific expectations."""
    k = interaction.parameters['k']
    eq = interaction.parameters['eq']

    # Force and energy should vanish at equilibrium
    F_eq = interaction.F(eq)
    U_eq = interaction.U(eq)
    assert isclose(F_eq, 0.0, abs_tol=1e-10), (
        f'{interaction.name}: F(eq) != 0, got {F_eq}'
    )
    assert isclose(U_eq, 0.0, abs_tol=1e-10), (
        f'{interaction.name}: U(eq) != 0, got {U_eq}'
    )

    # Symmetry around equilibrium
    dx = 0.05 * eq if eq > 0 else 0.05
    U_plus = interaction.U(eq + dx)
    U_minus = interaction.U(eq - dx)
    assert isclose(U_plus, U_minus, rel_tol=1e-10, abs_tol=1e-10), (
        f'{interaction.name}: U not symmetric about eq; '
        f'U(eq+dx)={U_plus}, U(eq-dx)={U_minus}'
    )

    # Restoring force direction
    F_right = interaction.F(eq + dx)
    F_left = interaction.F(eq - dx)
    assert F_right < 0.0, f'{interaction.name}: F(eq+dx) not restoring (<0), got {F_right}'
    assert F_left > 0.0, f'{interaction.name}: F(eq-dx) not restoring (>0), got {F_left}'

    # --- Water-specific numerical sanity checks ---
    # Distances/angles are in nm / radians, matching the example models.

    # O-H bond: eq ~ 0.1 nm (1 Å), energy near 0 there.
    if is_oh_bond(interaction):
        eq_oh = eq
        assert 0.09 <= eq_oh <= 0.11, (
            f"{interaction.name}: O-H equilibrium distance {eq_oh:.5f} nm outside "
            'expected water range [0.09, 0.11] nm'
        )

        x_ref = 0.10  # 1 Å in nm
        U_ref = interaction.U_base(x_ref)
        # Allow a few kJ/mol spread across force fields
        assert abs(U_ref) < 20.0, (
            f"{interaction.name}: bond energy at 0.10 nm should be near 0, got U={U_ref:.3f}"
        )

    # H-O-H angle: ~104.5° (~1.824 rad), energy near 0 there.
    if is_hoh_angle(interaction):
        eq_angle = eq
        assert 1.7 <= eq_angle <= 2.1, (
            f"{interaction.name}: H-O-H equilibrium angle {eq_angle:.4f} rad outside "
            'expected water range [1.7, 2.1] rad'
        )

        theta_ref = 1.824  # ~104.5 degrees in radians
        U_ref = interaction.U_base(theta_ref)
        assert abs(U_ref) < 10.0, (
            f"{interaction.name}: angle energy at ~104.5° should be near 0, got U={U_ref:.3f}"
        )


def sense_check_coulombic(interaction):
    """Sense-check for Coulombic interactions with water-like distances."""
    cp = interaction.parameters['charge_product']
    ecf = interaction.constants['ecf']
    epsilon_r = interaction.constants['epsilon_r']

    # If the net product is effectively zero, just rely on generic checks
    if isclose(cp, 0.0, abs_tol=1e-12):
        return

    cutoff = interaction.cutoff

    # Choose a typical intermolecular separation and a much shorter one
    if cutoff is not None:
        x_mid = min(0.3, 0.5 * cutoff)  # nm
        x_short = max(0.05, 0.1 * x_mid)
    else:
        x_mid = 0.3
        x_short = 0.05

    # Ensure we stay inside the cutoff for the base potential tests
    if cutoff is not None and x_mid >= cutoff:
        x_mid = 0.9 * cutoff
    if cutoff is not None and x_short >= cutoff:
        x_short = 0.5 * x_mid

    U_mid = interaction.U_base(x_mid)
    U_short = interaction.U_base(x_short)

    # Sign of the interaction at a typical distance
    if cp > 0:
        assert U_mid > 0.0, (
            f'{interaction.name}: like-charge U(x_mid) should be > 0, got {U_mid}'
        )
    else:
        assert U_mid < 0.0, (
            f'{interaction.name}: unlike-charge U(x_mid) should be < 0, got {U_mid}'
        )

    # "Very high" magnitude at short distances compared to a typical separation
    if not isclose(U_mid, 0.0, abs_tol=1e-12):
        assert abs(U_short) > 5.0 * abs(U_mid), (
            f'{interaction.name}: |U| at very short distance should be much higher '
            f'than at a typical distance; got U_short={U_short:.3f}, U_mid={U_mid:.3f}'
        )


def sense_check_lj(interaction):
    """Sense-check for Lennard-Jones interactions."""
    epsilon = interaction.parameters['epsilon']
    sigma = interaction.parameters['sigma']
    cutoff = interaction.cutoff

    # Use the shifted interaction (with cutoff) for the qualitative checks
    # Minimum at x = 2^(1/6) * sigma
    x_min = 2.0 ** (1.0 / 6.0) * sigma
    U_min = interaction.U(x_min)

    # Well should be negative near its minimum
    assert U_min < 0.0, (
        f'{interaction.name}: LJ U at the nominal minimum should be negative, got {U_min}'
    )

    # Repulsive at short range
    U_short = interaction.U(0.8 * sigma)
    assert U_short > 0.0, (
        f'{interaction.name}: LJ U at x<sigma should be positive (repulsive), got {U_short}'
    )

    # Attractive just beyond the minimum (inside cutoff if present)
    x_attr = 1.2 * x_min
    if cutoff is not None:
        x_attr = min(x_attr, 0.8 * cutoff)
    U_attr = interaction.U(x_attr)
    assert U_attr < 0.0, (
        f'{interaction.name}: LJ U just beyond minimum should be negative (attractive), got {U_attr}'
    )

    # As we move further out (but still inside cutoff), |U| should shrink,
    # approaching zero before the cutoff.
    x_far = 3.0 * sigma
    if cutoff is not None:
        x_far = min(x_far, 0.9 * cutoff)
        if x_far <= x_attr:
            x_far = 0.9 * cutoff

    U_far = interaction.U(x_far)
    assert abs(U_far) < abs(U_attr), (
        f'{interaction.name}: |U(x_far)| should be smaller than |U(x_attr)|; '
        f'got U_far={U_far}, U_attr={U_attr}'
    )


def sense_check_interaction(interaction, data):
    """Sense-check a single interaction using the physical criteria we agreed."""
    sense_check_common(interaction, data)

    potential = interaction.potential

    if potential is Harmonic:
        sense_check_harmonic(interaction)
    elif potential is Coulombic:
        sense_check_coulombic(interaction)
    elif potential is LennardJones:
        sense_check_lj(interaction)
    # For any other potentials we only apply the generic checks above.


def run_all():
    for ff_type in FORCE_FIELDS:
        ff = ff_type()
        ff_name = getattr(ff, 'name', ff_type.__name__)
        print(f'\n=== Testing interaction potentials in {ff_name} ===')

        assert ff.interaction_types, f'{ff_name} defines no interaction types'

        for interaction in ff.interaction_types:
            print(f'  - {interaction.name}')
            data = sample_interaction(interaction)
            sense_check_interaction(interaction, data)

    print('\nAll potential interaction tests passed.')


if __name__ == '__main__':
    run_all()
