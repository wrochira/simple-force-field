from dataclasses import dataclass
from abc import ABC, ABCMeta, abstractmethod

from .atom_type import AtomType


class PotentialMeta(ABCMeta):
    """Metaclass that enforces required fields and methods are implemented at class definition time."""
    def __new__(mcs, name, bases, namespace, **kwargs):
        cls = super().__new__(mcs, name, bases, namespace, **kwargs)
        # Check if this is a concrete class (not the base Potential class)
        if bases and bases[0].__name__ == 'Potential':
            # Check required fields and methods
            missing = [ ]
            for name in ('name', 'constant_names', 'parameter_names', 'F', 'U'):
                if name not in namespace:
                    missing.append(name)
            if missing:
                missing_str = ', '.join(missing)
                raise TypeError(f'Cannot create class "{name}": {missing_str} is required')
        return cls


@dataclass
class Potential(ABC, metaclass=PotentialMeta):
    """Base class for potential energy functions.

    All potentials must define:
    - name: str - The name of the potential
    - constant_names: frozenset[str] - The names of the constants of the potential
    - parameter_names: frozenset[str] - The names of the parameters of the potential
    - F() - Method to calculate force
    - U() - Method to calculate potential energy
    - calculate_parameters() - Method to calculate the parameters of the potential
    """
    name: str
    constant_names:  frozenset[str]
    parameter_names: frozenset[str]

    @abstractmethod
    def F(x: float, *args, **kwargs) -> float:
        """Calculate force (negative gradient of potential energy)."""
        raise NotImplementedError

    @abstractmethod
    def U(x: float, *args, **kwargs) -> float:
        """Calculate potential energy."""
        raise NotImplementedError

    @staticmethod
    def calculate_parameters(atom_types: tuple[AtomType, ...]) -> dict[str, float]:
        """Calculate the parameters of the potential."""
        raise NotImplementedError(f'This potential does not support parameter calculation')
