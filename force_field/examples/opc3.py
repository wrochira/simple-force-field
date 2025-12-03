"""
OPC3 water model

Parameter reference:
- J Chem Phys. 2016 Aug 15;145(7):074501. doi: 10.1063/1.4960175
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4991989/
"""

from ..core import ForceField, DecompositionType
from ..potentials import Harmonic, Coulombic, LennardJones


class OPC3(ForceField):
    def __init__(self):
        super().__init__()

        self.name = 'OPC3'
        self.flexible = False

        # Define atom types
        self.define_atom_type(
            name='Oxygen',
            symbol='OW',
            mass=15.9994,
            charge=-0.8952,
            sigma=0.31743,
            epsilon=0.68369
        )

        self.define_atom_type(
            name='Hydrogen',
            symbol='HW',
            mass=1.008,
            charge=+0.4476
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
