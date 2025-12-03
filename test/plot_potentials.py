"""
Generate visualisation plots for force field interaction potentials.

This script samples and visualises the force (F) and potential energy (U)
functions for all interactions across multiple water force field models.

Output structure:
- output/potentials/<force_field_name>/<interaction_name>.csv
  Individual CSV files with sampled r, F(r), U(r) data

- output/potentials/<force_field_name>/<interaction_name>.pdf
  Individual plots showing F(r) and U(r) for each interaction

- output/potentials/<force_field_name>/summary_potentials.pdf
  Summary plot showing all interactions within a single force field

- output/potentials/megasummary_potentials.pdf
  Mega-summary comparing the same interaction type across all force fields

The script automatically determines appropriate distance ranges for each
interaction based on parameters (equilibrium distance, sigma, cutoff) to
ensure physically relevant sampling ranges.

Usage:
    python3 plot_potentials.py

Requirements:
    matplotlib
"""

import sys
import csv
import math
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt

# Add parent directory to path to import force_field module
sys.path.insert(0, str(Path(__file__).parent.parent))

from force_field import (
    FBA,
    OPC3,
    SPC,
    SPC_E,
    SPC_FW,
    TIP3P,
    TIP3P_CHARMM,
    TIP3P_EW,
    TIP3P_FB,
    TIP3P_FS,
    TIP3P_PW,
)

FORCE_FIELDS = [
    FBA,
    OPC3,
    SPC,
    SPC_E,
    SPC_FW,
    TIP3P,
    TIP3P_CHARMM,
    TIP3P_EW,
    TIP3P_FB,
    TIP3P_FS,
    TIP3P_PW,
]

BASE_OUTPUT_DIR = Path(__file__).parent / 'output' / 'potentials'


def safe_name(name: str) -> str:
    """Make a filesystem-safe name from an interaction or force-field name."""
    return ''.join(c if c.isalnum() or c in ('_', '-', '.') else '_' for c in name)


def get_distance_range(interaction, default_min: float = 0.1, default_max: float = 1.0) -> tuple[float, float]:
    """Heuristically choose a sensible distance range for an interaction."""
    params = interaction.parameters
    cutoff = interaction.cutoff

    r_min = default_min
    r_max = default_max

    # Prefer explicit equilibrium distance if available (e.g. harmonic bonds)
    if 'eq' in params:
        eq = params['eq']
        r_min = 0.6 * eq
        r_max = 1.4 * eq

    # For Lennard-Jones style interactions, use sigma if available
    elif 'sigma' in params:
        sigma = params['sigma']
        r_min = 0.6 * sigma
        r_max = 2.5 * sigma

    # Otherwise, if a cutoff is defined, use a fraction of it
    elif cutoff is not None:
        r_min = 0.2 * cutoff
        r_max = 1.2 * cutoff

    # Ensure strictly positive and sane ordering
    r_min = max(r_min, 1e-3)
    if r_max <= r_min:
        r_max = r_min * 10.0

    return r_min, r_max


def sample_interaction(interaction, num_points: int = 200):
    """Sample F(r) and U(r) over a distance range for a single interaction."""
    r_min, r_max = get_distance_range(interaction)
    xs = [r_min + (r_max - r_min) * i / (num_points - 1) for i in range(num_points)]

    data = []
    for x in xs:
        F_val = interaction.F(x)
        U_val = interaction.U(x)
        data.append((x, F_val, U_val))

    return data


def write_csv(path: Path, data):
    """Write sampled interaction data to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['r', 'F', 'U'])
        writer.writerows(data)


def plot_interaction(path: Path, data, title: str):
    """Create a PDF with U(r) and F(r) for a single interaction."""
    xs = [d[0] for d in data]
    Fs = [d[1] for d in data]
    Us = [d[2] for d in data]

    fig, (ax_u, ax_f) = plt.subplots(2, 1, sharex=True, figsize=(6, 6))

    ax_u.plot(xs, Us, label='U(r)')
    ax_u.set_ylabel('U')
    ax_u.set_title(title)
    ax_u.grid(True, alpha=0.3)

    ax_f.plot(xs, Fs, label='F(r)', color='tab:red')
    ax_f.set_xlabel('r')
    ax_f.set_ylabel('F')
    ax_f.grid(True, alpha=0.3)

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def plot_summary(path: Path, interaction_results, title: str):
    """Create a summary PDF with separate subplots for each interaction in a force field."""
    n = len(interaction_results)
    if n == 0:
        return

    ncols = min(3, n)
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), sharex=False)
    if not isinstance(axes, (list, tuple)):
        # axes is a numpy array or a single Axes
        try:
            axes = axes.flatten()
        except AttributeError:
            axes = [axes]

    for ax, (name, data) in zip(axes, interaction_results):
        xs = [d[0] for d in data]
        Us = [d[2] for d in data]
        ax.plot(xs, Us)
        ax.set_title(name, fontsize='small')
        ax.set_xlabel('r')
        ax.set_ylabel('U')
        ax.grid(True, alpha=0.3)

    # Hide any unused subplots
    for ax in axes[n:]:
        ax.axis('off')

    fig.suptitle(title, fontsize='medium')
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def plot_megasummary(path: Path, interaction_groups):
    """Create a mega-summary PDF: one subplot per interaction name, overlaying force-field potentials."""
    interaction_names = sorted(interaction_groups.keys())
    n = len(interaction_names)
    if n == 0:
        return

    ncols = min(3, n)
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), sharex=False)
    if not isinstance(axes, (list, tuple)):
        try:
            axes = axes.flatten()
        except AttributeError:
            axes = [axes]

    for ax, name in zip(axes, interaction_names):
        for ff_name, data in interaction_groups[name]:
            xs = [d[0] for d in data]
            Us = [d[2] for d in data]
            ax.plot(xs, Us, label=ff_name)

        ax.set_title(name, fontsize='small')
        ax.set_xlabel('r')
        ax.set_ylabel('U')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize='x-small')

    # Hide any unused subplots
    for ax in axes[len(interaction_names):]:
        ax.axis('off')

    fig.suptitle('Interaction potentials across force fields', fontsize='medium')
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def main():
    BASE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Aggregate interaction data across all force fields for the mega-summary
    all_interaction_groups = defaultdict(list)

    for ff_type in FORCE_FIELDS:
        ff = ff_type()
        ff_name = getattr(ff, 'name', ff_type.__name__)
        ff_dir = BASE_OUTPUT_DIR / safe_name(ff_name)
        ff_dir.mkdir(parents=True, exist_ok=True)

        print(f'Generating data for force field: {ff_name}')

        interaction_results = []

        for interaction in ff.interaction_types:
            interaction_name = interaction.name
            print(f'  - Sampling interaction: {interaction_name}')

            data = sample_interaction(interaction)
            interaction_results.append((interaction_name, data))

            # Store for per-force-field outputs
            base = safe_name(interaction_name)
            csv_path = ff_dir / f'{base}.csv'
            pdf_path = ff_dir / f'{base}.pdf'

            write_csv(csv_path, data)
            plot_interaction(pdf_path, data, f'{ff_name} - {interaction_name}')

            # Store for mega-summary (grouped by interaction name)
            all_interaction_groups[interaction_name].append((ff_name, data))

        # Summary plot per force field (subplots per interaction)
        if interaction_results:
            summary_pdf = ff_dir / 'summary_potentials.pdf'
            plot_summary(summary_pdf, interaction_results, f'{ff_name} - interaction potentials')

    # Global mega-summary plot (subplots per interaction name, curves per force field)
    if all_interaction_groups:
        mega_pdf = BASE_OUTPUT_DIR / 'megasummary_potentials.pdf'
        plot_megasummary(mega_pdf, all_interaction_groups)


if __name__ == '__main__':
    main()
