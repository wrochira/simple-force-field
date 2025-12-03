from dataclasses import dataclass

from ..core import Potential


@dataclass
class Coulombic(Potential):
    name: str = 'coulombic'

    constant_names:  frozenset[str] = frozenset([ 'ecf', 'epsilon_r' ])
    parameter_names: frozenset[str] = frozenset([ 'charge_product' ])

    def F(x: float, charge_product: float, ecf: float, epsilon_r: float) -> float:
        return ecf * charge_product / (epsilon_r * x**2)

    def U(x: float, charge_product: float, ecf: float, epsilon_r: float) -> float:
        return ecf * charge_product / (epsilon_r * x)

    def calculate_parameters(atom_types):
        # Calculate the charge product of the two atom types
        charge_product = atom_types[0].charge * atom_types[1].charge
        return { 'charge_product': charge_product }
