from .atom_type import AtomType
from .potential import Potential
from ..constants import Constants


class InteractionType():
    """
    Defines an interaction between specific atom types with a potential function.

    An InteractionType represents a specific type of interaction (e.g., Lennard-Jones
    between oxygen atoms, Coulomb between charged particles) that can occur in a
    molecular system. It combines atom types, a potential function, force decomposition
    method, and the necessary parameters and constants.

    Attributes:
        name (str): Unique identifier for this interaction type.
        category_name (str): Category of the interaction
        bonded (bool): Whether this is a bonded interaction.
        potential (Potential): The potential function defining force and energy calculations.
        decomposition (str): Label for how forces should be decomposed into vector components.
        atom_types (tuple[AtomType, ...]): Tuple of atom types involved in this interaction.
        parameters (dict[str, float]): Parameters for the potential function (e.g., sigma, epsilon).
        constants (dict[str, float]): Physical constants required by the potential (e.g., ecf).
        cutoff (float, optional): Distance beyond which the interaction is set to zero.

    The class automatically:
    - Loads parameters from explicit values or calculates them from atom type properties
    - Loads constants from the Constants class or explicit values
    - Validates that all required parameters and constants are present
    """
    def __init__(self,
        name: str,
        category_name: str,
        bonded: bool,
        potential: Potential,
        decomposition: str,
        atom_types: tuple[AtomType, ...],
        parameters: dict[str, float] = None,
        constants: dict[str, float] = None,
        cutoff: float = None
    ):
        self.name = name
        self.category_name = category_name
        self.bonded = bonded
        self.potential = potential
        self.decomposition = decomposition
        self.atom_types = atom_types
        self.parameters = { }
        self.constants = { }
        self.cutoff = cutoff

        # Load parameters
        if parameters is None:
            parameters = { }
        for param_name in self.potential.parameter_names:
            if param_name in parameters:
                self.parameters[param_name] = parameters[param_name]
                continue
            try:
                self.parameters[param_name] = self.potential.calculate_parameters(self.atom_types)[param_name]
            except NotImplementedError:
                raise ValueError(f'Potential {self.potential.name} requires parameters be specified explicitly')
            except Exception as e:
                raise ValueError(f'Error calculating parameters for potential {self.potential.name}: {e}')
        assert all(param_name in self.parameters for param_name in self.potential.parameter_names), \
        f'Parameters {self.parameters} do not match potential {self.potential.name}'

        # Load constants
        if constants is None:
            constants = { }
        for const_name in self.potential.constant_names:
            if const_name in constants:
                self.constants[const_name] = constants[const_name]
                continue
            try:
                self.constants[const_name] = getattr(Constants, const_name)
            except AttributeError:
                raise ValueError(f'Constant {const_name} is not defined')
            except Exception as e:
                raise ValueError(f'Error getting value for constant {const_name}: {e}')
        assert all(const_name in self.constants for const_name in self.potential.constant_names), \
        f'Constants {self.constants} do not match potential {self.potential.name}'


    def F_base(self, x: float) -> float:
        """Calculate force without cutoff handling."""
        return self.potential.F(x, **self.parameters, **self.constants)


    def U_base(self, x: float) -> float:
        """Calculate potential energy without cutoff handling."""
        return self.potential.U(x, **self.parameters, **self.constants)


    def F(self, x: float) -> float:
        """Calculate force with cutoff handling."""
        if self.cutoff is not None and x > self.cutoff:
            return 0.0
        return self.F_base(x)


    def U(self, x: float) -> float:
        """Calculate potential energy with cutoff handling."""
        if self.cutoff is not None:
            if x > self.cutoff:
                return 0.0
            return self.U_base(x) - self.U_base(self.cutoff)
        return self.U_base(x)
