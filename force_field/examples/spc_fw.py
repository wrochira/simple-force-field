"""
SPC/Fw water model

Parameter reference:
- TODO: add reference
"""

from ..core import ForceField, DecompositionType
from ..potentials import Harmonic, Coulombic, LennardJones


class SPC_FW(ForceField):
    def __init__(self):
        super().__init__()

        self.name = 'SPC/Fw'
        self.flexible = True

        # Define atom types
        self.define_atom_type(
            name='Oxygen',
            symbol='OW',
            mass=15.9994,
            charge=-0.82,
            sigma=0.31655,
            epsilon=0.6503
        )

        self.define_atom_type(
            name='Hydrogen',
            symbol='HW',
            mass=1.008,
            charge=+0.41
        )

        # Define interaction types
        self.define_interaction_type(
            name='O-H Bond',
            category_name='Bond Stretching',
            potential=Harmonic,
            decomposition=DecompositionType.bond,
            parameters={'k': 423_664, 'eq': 0.1012},
            bonded=True,
            atom_types=('OW', 'HW')
        )

        self.define_interaction_type(
            name='H-O-H Angle',
            category_name='Angle Bending',
            potential=Harmonic,
            decomposition=DecompositionType.angle,
            parameters={'k': 317.56, 'eq': 1.97641},
            bonded=True,
            atom_types=('HW', 'OW', 'HW')
        )

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
