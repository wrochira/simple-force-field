<p align="center">
    <img src=".assets/img/icon.svg" alt="Project Icon" width="150">
    <br>
    <h1 align="center">
        Simple Force Field
    </h1>
    <h3 align="center">
        A simple, extensible library for defining and working with molecular dynamics force fields in Python.
    </h3>
    <p></p>
    <p align="center">
    <a href="#about">About</a> •
    <a href="#key-features">Key Features</a> •
    <a href="#quick-start">Quick Start</a> •
    <a href="#learn-more">Learn More</a>
    </p>
    <p></p>
</p>

## About

`force_field` is a Python library for defining molecular dynamics force fields from scratch. Rather than wrapping complex simulation engines or relying on external parameter files, it provides a clean, minimal framework for building force field definitions programmatically. You define atom types, potentials, and interaction rules using pure Python, and the library handles parameter derivation, cutoff behaviour, and force/energy calculations.

The library was designed for clarity and extensibility. It's particularly well‑suited for:

- Learning how force fields work by implementing them yourself
- Rapid prototyping of new potentials or water models
- Generating force field definitions for custom MD engines
- Educational projects requiring transparent, inspectable force field logic

All currently implemented "force fields" are three‑site water models (SPC, SPC/E, SPC/Fw, TIP3P, and several TIP3P variants). The architecture is general enough to support macromolecular force fields, torsional potentials, and bonded interactions, but these are not yet implemented.

## Key Features

- **Clean, modular architecture** - Atom types, potentials, and interaction types are separate, reusable components
- **Extensible potential system** - Define new potentials by implementing force (`F`) and energy (`U`) methods; the library handles the rest
- **Automatic parameter derivation** - Parameters like Lennard‑Jones $\sigma$ and $\epsilon$ can be calculated automatically from atom type properties using combination rules
- **Cutoff support** - Interactions can specify a cutoff distance; forces and energies are automatically shifted to zero beyond the cutoff
- **No special treatment for water** - Water models are implemented as standard force fields using the same abstractions as any other molecular system
- **Pre‑defined water models** - Includes eleven three‑site water models: FBA, OPC3, SPC, SPC/E, SPC/Fw, TIP3P, TIP3P‑CHARMM, TIP3P‑EW, TIP3P‑FB, TIP3P‑FS, and TIP3P‑PW
- **Compatible with Python MD libraries** - Designed to work seamlessly with other Python MD tools (but not the standalone C MD engine modules)
- **Well‑tested** - Includes comprehensive test suites covering basic structure, potential values, and cutoff behaviour

## Quick Start

Install the dependencies:

```bash
pip install numpy pytest matplotlib
```

Use a pre‑defined water model:

```python
from force_field import TIP3P

# Instantiate the force field
ff = TIP3P()

# Access atom types
oxygen = ff.get_atom_type('OW')
print(f"Oxygen mass: {oxygen.mass} AMU")
print(f"Oxygen charge: {oxygen.charge} e")

# Access interaction types and calculate forces/energies
for interaction in ff.interaction_types:
    print(f"{interaction.name}: F(1.0 nm) = {interaction.F(1.0):.3f} kJ/(mol·nm)")
```

Define your own simple force field:

```python
from force_field import ForceField, DecompositionType
from force_field.potentials import LennardJones

ff = ForceField(name='My Force Field')

# Define atom types
ff.define_atom_type(name='Carbon', symbol='C', mass=12.011, sigma=0.34, epsilon=0.36)
ff.define_atom_type(name='Nitrogen', symbol='N', mass=14.007, sigma=0.325, epsilon=0.71)

# Define interactions
ff.define_interaction_type(
    name='C-N Dispersion',
    category_name='Dispersion',
    bonded=False,
    potential=LennardJones,
    decomposition=DecompositionType.pairwise,
    atom_types=('C', 'N'),
    cutoff=1.0
)

# The library automatically calculates σ and ε for the C-N pair
# using the Lorentz-Berthelot combination rules
```

Run the tests:

```bash
pytest
```

## Learn More

For detailed documentation, including:

- System architecture and design philosophy
- How to define new potentials with custom force and energy functions
- How to define new force fields and water models
- Complete API reference for all core classes
- Example implementations and usage patterns
- Testing guidance

See **[DOCUMENTATION.md](DOCUMENTATION.md)**.

## Status

This library is functional and tested, but not under active development aside from occasional refinements.

## TODO

- Ability to combine force fields (e.g. a macromolecular force field and a water model)
- Implement macromolecular force fields and torsional potentials
- Finish the SPC/E model by adding an average polarisation correction to the potential energy function
    * $ E_{\text{pol}}={\frac {1}{2}}\sum _{i}{\frac {(\mu -\mu ^{0})^{2}}{\alpha _{i}}} $
- Complete citations, sources, and descriptions in the force field files

## License

This project is licensed under the [CC BY-NC-SA 4.0 Licence](https://creativecommons.org/licenses/by-nc-sa/4.0/).
You can freely use, modify, and share this work for non-commercial purposes. You must give appropriate credit.
Any modifications must be shared under the same licence. You cannot use this work for commercial purposes without permission.

[![CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
