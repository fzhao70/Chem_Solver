#!/usr/bin/env python3
"""
Example: Single Grid Cell

Demonstrates solving chemistry for a single grid cell (same as multi-step example).
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from chem_solver import ChemicalSolver


def main():
    print("=" * 60)
    print("SMVGEAR-like Solver: Single Grid Cell Example")
    print("=" * 60)

    # Initialize solver
    mechanism_file = os.path.join(os.path.dirname(__file__), '..', 'mechanisms', 'simple_chemistry.json')
    solver = ChemicalSolver(mechanism_file, time_step=60.0)

    # Single grid cell initial conditions
    concentrations = {
        'O3': 1.0e12,
        'NO': 1.0e11,
        'NO2': 5.0e11,
        'O1D': 0.0,
        'O3P': 0.0,
        'M': 2.5e19
    }

    temperature = 298.0
    pressure = 101325.0
    photolysis_rates = {
        'JNO2': 0.008,
        'JO3_O1D': 2.0e-5
    }

    total_time = 3600.0  # 1 hour

    print("\nSingle Grid Cell Simulation")
    print(f"  Integration time: {total_time / 60.0} minutes")
    print(f"  Time step: {solver.time_step} s")

    print("\nInitial O3: {:.2e} molecules/cm³".format(concentrations['O3']))
    print("Initial NO: {:.2e} molecules/cm³".format(concentrations['NO']))
    print("Initial NO2: {:.2e} molecules/cm³".format(concentrations['NO2']))

    # Solve for single grid cell
    final = solver.solve(
        concentrations,
        temperature,
        pressure,
        photolysis_rates,
        total_time=total_time
    )

    print("\nFinal O3: {:.2e} molecules/cm³".format(final['O3']))
    print("Final NO: {:.2e} molecules/cm³".format(final['NO']))
    print("Final NO2: {:.2e} molecules/cm³".format(final['NO2']))

    print("\n" + "=" * 60)


if __name__ == '__main__':
    main()
