from dataclasses import dataclass

from ..core import Potential


@dataclass
class LennardJones(Potential):
    name: str = 'lennard_jones'

    constant_names:  frozenset[str] = frozenset()
    parameter_names: frozenset[str] = frozenset([ 'epsilon', 'sigma' ])

    def F(x: float, epsilon: float, sigma: float) -> float:
        return 48 * epsilon * ((sigma**12 / x**13) - 0.5 * (sigma**6 / x**7))

    def U(x: float, epsilon: float, sigma: float) -> float:
        return 4 * epsilon * ((sigma / x)**12 - (sigma / x)**6)

    def calculate_parameters(atom_types):
        # Calculate epsilon using the Lorentz-Berthelot combination rule (geometric mean)
        epsilon = (atom_types[0].epsilon * atom_types[1].epsilon) ** 0.5
        # Calculate sigma using the Lorentz-Berthelot combination rule (arithmetic mean)
        sigma = (atom_types[0].sigma + atom_types[1].sigma) / 2
        return { 'epsilon': epsilon, 'sigma': sigma }
