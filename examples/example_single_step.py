#!/usr/bin/env python3
"""
Example: Single Step Solver

Demonstrates the use of the one-step solver API.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from chem_solver import ChemicalSolver


def main():
    print("=" * 60)
    print("SMVGEAR-like Solver: Single Step Example")
    print("=" * 60)

    # Initialize solver with simple tropospheric chemistry
    mechanism_file = os.path.join(os.path.dirname(__file__), '..', 'mechanisms', 'simple_chemistry.json')
    solver = ChemicalSolver(mechanism_file, time_step=60.0)  # 60 second time step

    # Initial concentrations (molecules/cm^3)
    initial_concentrations = {
        'O3': 1.0e12,      # 1e12 molecules/cm^3 (~40 ppb at STP)
        'NO': 1.0e11,      # 1e11 molecules/cm^3 (~4 ppb at STP)
        'NO2': 5.0e11,     # 5e11 molecules/cm^3 (~20 ppb at STP)
        'O1D': 0.0,
        'O3P': 0.0,
        'M': 2.5e19        # Air density at STP
    }

    # Environmental conditions
    temperature = 298.0  # K (25°C)
    pressure = 101325.0  # Pa (1 atm)

    # Photolysis rates (1/s)
    photolysis_rates = {
        'JNO2': 0.008,      # NO2 photolysis rate
        'JO3_O1D': 2.0e-5   # O3 -> O1D photolysis rate
    }

    print("\nInitial Conditions:")
    print(f"  Temperature: {temperature} K")
    print(f"  Pressure: {pressure} Pa")
    print(f"  Time step: {solver.time_step} s")
    print("\nInitial Concentrations:")
    for species, conc in initial_concentrations.items():
        if species != 'M':
            print(f"  {species:8s}: {conc:.2e} molecules/cm³")

    # Perform single step
    print("\nPerforming single integration step...")
    final_concentrations = solver.solve_step(
        initial_concentrations,
        temperature,
        pressure,
        photolysis_rates
    )

    print("\nFinal Concentrations (after one step):")
    for species, conc in final_concentrations.items():
        if species != 'M':
            initial = initial_concentrations.get(species, 0)
            change = ((conc - initial) / max(initial, 1e-30)) * 100 if initial > 0 else 0
            print(f"  {species:8s}: {conc:.2e} molecules/cm³ (change: {change:+.2f}%)")

    print("\n" + "=" * 60)


if __name__ == '__main__':
    main()
