from dataclasses import dataclass

from ..core import Potential


@dataclass
class Harmonic(Potential):
    name: str = 'harmonic'

    constant_names:  frozenset[str] = frozenset()
    parameter_names: frozenset[str] = frozenset([ 'k', 'eq' ])

    def F(x: float, k: float, eq: float) -> float:
        return -k * (x - eq)

    def U(x: float, k: float, eq: float) -> float:
        return 0.5 * k * (x - eq) ** 2
