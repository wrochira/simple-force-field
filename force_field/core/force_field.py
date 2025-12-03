from itertools import combinations

from .atom_type import AtomType
from .potential import Potential
from .decomposition import DecompositionType
from .interaction_type import InteractionType


class ForceField:
    """
    A force field definition for molecular dynamics simulations.

    The ForceField class manages atom types and interaction types that define
    the forces between atoms in a molecular system. It provides methods to
    define atom types with their properties and interaction types with their
    potentials, parameters, and decomposition methods.

    Attributes:
        name (str): The name of the force field.
        atom_types (list[AtomType]): List of defined atom types.
        atom_type_identifiers (set[str]): Set of all atom type names and symbols
            for quick lookup and uniqueness checking.
        interaction_types (list[InteractionType]): List of defined interaction types.
        interaction_type_identifiers (set[str]): Set of all interaction type names
            for quick lookup and uniqueness checking.
    """
    def __init__(self, name: str = None):
        self.name = name
        self.atom_types = [ ]
        self.atom_type_identifiers = set()
        self.interaction_types = [ ]
        self.interaction_type_identifiers = set()


    def define_atom_type(
        self,
        name: str,
        symbol: str,
        mass: float,
        **kwargs
    ):
        """Define an atom type."""
        assert name != symbol, 'Name and symbol cannot be the same'
        assert name not in self.atom_type_identifiers, \
            f'AtomType name {name} already defined as a name or symbol'
        assert symbol not in self.atom_type_identifiers, \
            f'AtomType symbol {symbol} already defined as a name or symbol'
        atom_type = AtomType(name=name, symbol=symbol, mass=mass, **kwargs)
        atom_type.interaction_types = [ ]
        self.atom_types.append(atom_type)
        self.atom_type_identifiers.add(name)
        self.atom_type_identifiers.add(symbol)


    def get_atom_type(self, identifier: str) -> AtomType:
        """Get an atom type by name or symbol."""
        if identifier in self.atom_type_identifiers:
            return next((atom_type for atom_type in self.atom_types if atom_type.name == identifier or atom_type.symbol == identifier), None)
        else:
            raise ValueError(f'AtomType type {identifier} not found')


    def define_interaction_type(
        self,
        name: str,
        category_name: str,
        bonded: bool,
        potential: Potential,
        decomposition: DecompositionType,
        atom_types: tuple[AtomType, str, ...],
        parameters: dict[str, float] = None,
        constants: dict[str, float] = None,
        cutoff: float = None
    ):
        """Define an interaction type."""

        assert name not in self.interaction_type_identifiers, f'Interaction type name {name} already defined'

        # Get atom type objects from identifiers
        atom_types_objects = [ ]
        try:
            for atom_type in atom_types:
                if isinstance(atom_type, AtomType):
                    atom_types_objects.append(atom_type)
                else:
                    atom_types_objects.append(self.get_atom_type(atom_type))
        except ValueError as e:
            raise ValueError(f'Atom type {atom_type} not found')
        except Exception as e:
            raise ValueError(f'Error getting atom types: {e}')
        atom_types = atom_types_objects

        # Create interaction type
        interaction_type = InteractionType(
            name=name,
            category_name=category_name,
            bonded=bonded,
            potential=potential,
            decomposition=decomposition,
            atom_types=atom_types,
            parameters=parameters,
            constants=constants,
            cutoff=cutoff
        )
        self.interaction_types.append(interaction_type)
        self.interaction_type_identifiers.add(name)
        for atom_type in set(atom_types):
            atom_type.interaction_types.append(interaction_type)
