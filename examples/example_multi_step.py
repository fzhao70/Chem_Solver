#!/usr/bin/env python3
"""
Example: Multi-Step Solver

Demonstrates the use of the multi-step solver API over a longer time period.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from chem_solver import ChemicalSolver


def main():
    print("=" * 60)
    print("SMVGEAR-like Solver: Multi-Step Example")
    print("=" * 60)

    # Initialize solver with methane oxidation chemistry
    mechanism_file = os.path.join(os.path.dirname(__file__), '..', 'mechanisms', 'methane_oxidation.json')
    solver = ChemicalSolver(mechanism_file, time_step=300.0)  # 5 minute time step

    # Initial concentrations (molecules/cm^3)
    initial_concentrations = {
        'CH4': 4.5e13,     # ~1.8 ppm at STP
        'OH': 1.0e6,       # Typical daytime OH
        'CH3O2': 0.0,
        'CH3OOH': 0.0,
        'HCHO': 0.0,
        'HO2': 1.0e8,
        'H2O': 5.0e17,
        'O2': 5.0e18,
        'M': 2.5e19
    }

    # Environmental conditions
    temperature = 298.0  # K
    pressure = 101325.0  # Pa

    # Photolysis rates (1/s)
    photolysis_rates = {
        'JCH3OOH': 5.0e-6,
        'JHCHO': 3.0e-5
    }

    # Simulation parameters
    total_time = 3600.0 * 24.0  # 24 hours
    time_step = 300.0           # 5 minutes

    print("\nInitial Conditions:")
    print(f"  Temperature: {temperature} K")
    print(f"  Pressure: {pressure} Pa")
    print(f"  Total time: {total_time / 3600.0} hours")
    print(f"  Time step: {time_step} s")

    print("\nInitial Concentrations:")
    for species, conc in initial_concentrations.items():
        if species not in ['M', 'O2', 'H2O']:
            print(f"  {species:8s}: {conc:.2e} molecules/cm³")

    # Solve over multiple steps
    print("\nSolving chemistry over 24 hours...")
    final_concentrations = solver.solve(
        initial_concentrations,
        temperature,
        pressure,
        photolysis_rates,
        total_time=total_time,
        time_step=time_step
    )

    print("\nFinal Concentrations (after 24 hours):")
    for species, conc in final_concentrations.items():
        if species not in ['M', 'O2', 'H2O']:
            initial = initial_concentrations.get(species, 0)
            change = ((conc - initial) / max(initial, 1e-30)) * 100 if initial > 0 else 0
            print(f"  {species:8s}: {conc:.2e} molecules/cm³ (change: {change:+.2f}%)")

    print("\n" + "=" * 60)


if __name__ == '__main__':
    main()
