"""
TIP3P-FB water model

Parameter reference:
- TODO: add reference
"""

from ..core import ForceField, DecompositionType
from ..potentials import Harmonic, Coulombic, LennardJones


class TIP3P_FB(ForceField):
    def __init__(self):
        super().__init__()

        self.name = 'TIP3P-FB'
        self.flexible = False

        # Define atom types
        self.define_atom_type(
            name='Oxygen',
            symbol='OW',
            mass=15.9994,
            charge=-0.8484,
            sigma=0.3178,
            epsilon=0.65214
        )

        self.define_atom_type(
            name='Hydrogen',
            symbol='HW',
            mass=1.008,
            charge=+0.4242
        )

        # Define interaction types
        self.define_interaction_type(
            name='H-H Electrostatic',
            category_name='Electrostatic',
            potential=Coulombic,
            decomposition=DecompositionType.pairwise,
            bonded=False,
            atom_types=('HW', 'HW'),
            cutoff=1.0
        )

        self.define_interaction_type(
            name='O-H Electrostatic',
            category_name='Electrostatic',
            potential=Coulombic,
            decomposition=DecompositionType.pairwise,
            bonded=False,
            atom_types=('OW', 'HW'),
            cutoff=1.0
        )

        self.define_interaction_type(
            name='O-O Electrostatic',
            category_name='Electrostatic',
            potential=Coulombic,
            decomposition=DecompositionType.pairwise,
            bonded=False,
            atom_types=('OW', 'OW'),
            cutoff=1.0
        )

        self.define_interaction_type(
            name='O-O Dispersion',
            category_name='Dispersion',
            potential=LennardJones,
            decomposition=DecompositionType.pairwise,
            bonded=False,
            atom_types=('OW', 'OW'),
            cutoff=1.0
        )
