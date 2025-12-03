from dataclasses import dataclass

@dataclass
class Constants:
    ecf: float        = 138.935458     # Electric conversion factor (kJ mol-1 nm e-2)
    epsilon_r: float  = 1.0            # Dielectric constant (unitless)
