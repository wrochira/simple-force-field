"""
TIP3P (pairwise) water model

Parameter reference:
- Rochira W. 2024
- These parameters were calculated using an ad-hoc parameter search algorithm
  that minimised errors in geometric and thermodynamic properties compared to the
  default pairwise TIP3P model. The goal was to create a model that emulates the
  behaviour of the standard TIP3P model without implementing an explicit angle-
  bending interaction - instead applying an H-H 'pseudobond'.
"""

from ..core import ForceField, DecompositionType
from ..potentials import Harmonic, Coulombic, LennardJones


class TIP3P_PW(ForceField):
    def __init__(self):
        super().__init__()

        self.name = 'TIP3P_PW'
        self.flexible = True

        # Define atom types
        self.define_atom_type(
            name='Oxygen',
            symbol='OW',
            mass=15.9994,
            charge=-0.834,
            sigma=0.31507,
            epsilon=0.63639
        )

        self.define_atom_type(
            name='Hydrogen',
            symbol='HW',
            mass=1.008,
            charge=+0.417,
        )

        # Define interaction types
        self.define_interaction_type(
            name='O-H Bond',
            category_name='Bond Stretching',
            potential=Harmonic,
            decomposition=DecompositionType.bond,
            parameters={'k': 502_416, 'eq': 0.09572},
            bonded=True,
            atom_types=('OW', 'HW')
        )

        self.define_interaction_type(
            name='H-H Pseudobond',
            category_name='Bond Stretching',
            potential=Harmonic,
            decomposition=DecompositionType.bond,
            parameters={'k': 145_330, 'eq': 0.15139},
            bonded=True,
            atom_types=('HW', 'HW')
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
