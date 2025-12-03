from dataclasses import dataclass


@dataclass(frozen=True)
class DecompositionType:
    """Decomposition type labels used by the force-field definitions.

    These are purely declarative tags; the molecular dynamics engine is
    responsible for interpreting them and performing the actual force
    decomposition in Cartesian coordinates.
    """

    pairwise: str = 'pairwise'
    bond: str = 'bond'
    angle: str = 'angle'
    dihedral: str = 'dihedral'
