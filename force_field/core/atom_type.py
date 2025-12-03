class AtomType:
    """Represents an atom type with required and optional attributes.

    Required attributes:
        name: str - Name of the atom type (e.g., 'O', 'H')
        symbol: str - Chemical symbol (e.g., 'O', 'H')
        mass: float - Atomic mass (AMU)

    Optional attributes can be passed as keyword arguments, e.g.
        charge: float - Atomic charge (e)
        sigma: float - Equilibrium distance for Lennard-Jones (nm)
        epsilon: float - Depth of potential well (kJ mol^-1)
        ...and any others
    """

    def __init__(self, name: str, symbol: str, mass: float, **kwargs):
        self.name = name
        self.symbol = symbol
        self.mass = mass
        self.interaction_types = None

        # Store all extra kwargs as attributes
        for key, value in kwargs.items():
            setattr(self, key, value)
